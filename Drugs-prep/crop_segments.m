function OUT = crop_segments(meta_file, scores_dir, mats_dir, out_dir, varargin)
% crop_baseline_ambtemp_segments (generalised)
% -------------------------------------------------------------------------
% For each row in meta_file (CSV or XLSX) with columns:
%   "Mouse" | "Genotype" | "Time started (s)" | "Time started (h)" [|
%   "Time finished (s)" | "Time finished (h)" (optional)]
%
% For that mouse & genotype, it:
%   - finds recordings for each requested condition:
%       YYYYMMDD_<cond>_mouseX_GENO_scored_scores_1Hz.csv
%       YYYYMMDD_<cond>_mouseX_GENO.mat
%   - if "Time finished (s)" exists:
%         crops them to [Time started (s), Time finished (s)]
%     else:
%         crops them to [Time started (s), END OF RECORDING]
%   - writes *_crop.csv and *_crop.mat into out_dir.
%
% Name–value options:
%   'eegVar'    : name of EEG vector in MAT      (default 'eeg')
%   'fsVar'     : name of sampling rate in MAT   (default 'eeg_frequency')
%   'ambVar'    : name of ambient temp vector    (default '', i.e. none)
%   'ambFsVar'  : sampling rate var for ambVar   (default '')
%   'conditions': string / cellstr / string array
%                 (default ["baseline","ambtemp"])
% -------------------------------------------------------------------------

p = inputParser;
addParameter(p,'eegVar','eeg',@(s)ischar(s)||isstring(s));
addParameter(p,'fsVar','eeg_frequency',@(s)ischar(s)||isstring(s));
addParameter(p,'ambVar','',@(s)ischar(s)||isstring(s));
addParameter(p,'ambFsVar','',@(s)ischar(s)||isstring(s));
addParameter(p,'conditions',["baseline","ambtemp"], ...
    @(c) isstring(c) || iscellstr(c) || ischar(c));

parse(p, varargin{:});

eegVar   = string(p.Results.eegVar);
fsVar    = string(p.Results.fsVar);
ambVar   = string(p.Results.ambVar);
ambFsVar = string(p.Results.ambFsVar);

conds_in = p.Results.conditions;
if ischar(conds_in)
    conds = string({conds_in});
elseif iscellstr(conds_in)
    conds = string(conds_in);
else
    conds = string(conds_in);
end

if nargin < 4 || isempty(out_dir)
    out_dir = fullfile(pwd,'cropped_segments');
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% ---------- Load meta table (CSV or XLSX; preserve original headers) -----
META = readtable(meta_file, 'VariableNamingRule','preserve');

reqCols = ["Mouse","Genotype","Time started (s)"];
if ~all(ismember(reqCols, META.Properties.VariableNames))
    error('Meta file must contain columns: %s', strjoin(reqCols, ', '));
end

MouseNum = META.("Mouse");
Genotype = string(META.("Genotype"));
t_start  = META.("Time started (s)");

% Optional "Time finished (s)" column:
if ismember("Time finished (s)", META.Properties.VariableNames)
    t_end = META.("Time finished (s)");
else
    t_end = nan(height(META),1);   % indicates "to end of recording"
end

OUT.entries = struct([]);
kEntry = 0;

for i = 1:numel(MouseNum)
    mouse_num   = MouseNum(i);
    mouse_label = "mouse" + string(mouse_num);
    geno        = Genotype(i);
    t0          = t_start(i);

    if isnan(t0)
        warning('Row %d: invalid start time, skipping.', i);
        continue;
    end

    % t1 may be NaN (meaning "to end")
    if i <= numel(t_end)
        t1 = t_end(i);
    else
        t1 = NaN;
    end

    has_t1 = ~isnan(t1);

    if has_t1 && t1 <= t0
        warning('Row %d: finish time <= start time, skipping.', i);
        continue;
    end

    if has_t1
        fprintf('\n=== Cropping %s (%s) from %.1f to %.1f s ===\n', ...
            mouse_label, geno, t0, t1);
    else
        fprintf('\n=== Cropping %s (%s) from %.1f s to END ===\n', ...
            mouse_label, geno, t0);
    end

    % For each requested condition
    for cond = conds

        fprintf('  Condition: %s\n', cond);

        % ---------- 1) Find scores CSV ----------
        scores_pattern = sprintf('*_%s_%s_%s_scored_scores_1Hz.csv', ...
                                 cond, mouse_label, geno);
        scores_file = find_one_file(scores_dir, scores_pattern);
        if scores_file == ""
            warning('    No scores CSV found matching "%s".', scores_pattern);
            continue;
        end

        % ---------- 2) Find MAT file ----------
        mat_pattern = sprintf('*_%s_%s_%s.mat', cond, mouse_label, geno);
        mat_file = find_one_file(mats_dir, mat_pattern);
        if mat_file == ""
            warning('    No MAT file found matching "%s".', mat_pattern);
            continue;
        end

        % ---------- 3) Crop scores CSV ----------
        T = readtable(scores_file);
        if ~ismember('time_s', T.Properties.VariableNames)
            error('Scores file %s has no "time_s" column.', scores_file);
        end

        if has_t1
            keep = (T.time_s >= t0) & (T.time_s <= t1);
        else
            keep = (T.time_s >= t0);     % to end of recording
        end

        T_crop = T(keep, :);

        % Rebase time so that manipulation starts at 0
        T_crop.time_s = T_crop.time_s - t0;

        [~, scores_base, ~] = fileparts(scores_file);
        out_scores_file = fullfile(out_dir, scores_base + "_crop.csv");
        writetable(T_crop, out_scores_file);
        fprintf('    -> Cropped scores: %s (n=%d rows)\n', ...
                out_scores_file, height(T_crop));

        % ---------- 4) Crop MAT signals ----------
        S = load(mat_file);

        if ~isfield(S, eegVar) || ~isfield(S, fsVar)
            warning('    MAT file %s missing %s or %s. Skipping MAT crop.', ...
                    mat_file, eegVar, fsVar);
            out_mat_file = "";
        else
            eeg    = S.(eegVar);
            fs_eeg = S.(fsVar);

            nSamp = numel(eeg);
            i0 = max(1, floor(t0 * fs_eeg) + 1);

            if has_t1
                i1 = min(nSamp, ceil(t1 * fs_eeg));
            else
                i1 = nSamp;   % to end of recording
            end

            eeg_crop = eeg(i0:i1);

            S_crop = S;         % start with original struct
            S_crop.(eegVar) = eeg_crop;
            S_crop.crop_t0_s = t0;

            % Actual absolute end time of the cropped segment
            crop_t1_eff = (i1 - 1) / fs_eeg;
            S_crop.crop_t1_s = crop_t1_eff;

            S_crop.crop_cond     = cond;
            S_crop.crop_mouse    = mouse_label;
            S_crop.crop_genotype = geno;

            % Optional ambient temperature channel
            if ambVar ~= "" && isfield(S, ambVar)
                amb = S.(ambVar);
                if ambFsVar ~= "" && isfield(S, ambFsVar)
                    fs_amb = S.(ambFsVar);
                else
                    fs_amb = fs_eeg; % assume same
                end
                nAmb = numel(amb);
                j0 = max(1, floor(t0 * fs_amb) + 1);

                if has_t1
                    j1 = min(nAmb, ceil(t1 * fs_amb));
                else
                    j1 = nAmb;
                end
                S_crop.(ambVar) = amb(j0:j1);
            end

            [~, mat_base, ~] = fileparts(mat_file);
            out_mat_file = fullfile(out_dir, mat_base + "_crop.mat");
            save(out_mat_file, '-struct', 'S_crop');
            fprintf('    -> Cropped MAT:    %s\n', out_mat_file);
        end

        % ---------- 5) Record entry ----------
        kEntry = kEntry + 1;
        OUT.entries(kEntry).mouse_num   = mouse_num;
        OUT.entries(kEntry).mouse_label = mouse_label;
        OUT.entries(kEntry).genotype    = geno;
        OUT.entries(kEntry).condition   = cond;
        OUT.entries(kEntry).t_start_s   = t0;
        OUT.entries(kEntry).t_end_s     = S_crop.crop_t1_s; % actual end
        OUT.entries(kEntry).scores_in   = scores_file;
        OUT.entries(kEntry).scores_out  = out_scores_file;
        OUT.entries(kEntry).mat_in      = mat_file;
        OUT.entries(kEntry).mat_out     = out_mat_file;
    end
end

fprintf('\nDone. Cropped segments are in: %s\n', out_dir);
end


% =====================================================================
function fullpath = find_one_file(dirpath, pattern)
% Helper: find exactly one file matching pattern in dirpath.
D = dir(fullfile(dirpath, pattern));
if isempty(D)
    fullpath = "";
elseif numel(D) > 1
    warning('    Multiple files match "%s" in %s; using the first (%s).', ...
            pattern, dirpath, D(1).name);
    fullpath = fullfile(D(1).folder, D(1).name);
else
    fullpath = fullfile(D(1).folder, D(1).name);
end
end
