function OUT = crop_baseline_ambtemp_segments(meta_csv, scores_dir, mats_dir, out_dir, varargin)
% crop_baseline_ambtemp_segments
% -------------------------------------------------------------------------
% For each row in meta_csv (ambtemp experiment):
%   Mouse | Genotype | Time started (s) | Time started (h) |
%                          Time finished (s) | Time finished (h)
%
% For that mouse & genotype, it:
%   - finds BASELINE and AMBTEMP recordings:
%       YYYYMMDD_baseline_mouseX_GENO_scored_scores_1Hz.csv
%       YYYYMMDD_ambtemp_mouseX_GENO_scored_scores_1Hz.csv
%       YYYYMMDD_baseline_mouseX_GENO.mat
%       YYYYMMDD_ambtemp_mouseX_GENO.mat
%   - crops them to [Time started (s), Time finished (s)]
%   - writes *_crop.csv and *_crop.mat into out_dir.
%
% Name–value options:
%   'eegVar'   : name of EEG vector in MAT      (default 'eeg')
%   'fsVar'    : name of sampling rate in MAT   (default 'eeg_frequency')
%   'ambVar'   : name of ambient temp vector    (default '', i.e. none)
%   'ambFsVar' : sampling rate var for ambVar   (default '')
% -------------------------------------------------------------------------

p = inputParser;
addParameter(p,'eegVar','eeg',@(s)ischar(s)||isstring(s));
addParameter(p,'fsVar','eeg_frequency',@(s)ischar(s)||isstring(s));
addParameter(p,'ambVar','',@(s)ischar(s)||isstring(s));
addParameter(p,'ambFsVar','',@(s)ischar(s)||isstring(s));
parse(p, varargin{:});

eegVar   = string(p.Results.eegVar);
fsVar    = string(p.Results.fsVar);
ambVar   = string(p.Results.ambVar);
ambFsVar = string(p.Results.ambFsVar);

if nargin < 4 || isempty(out_dir)
    out_dir = fullfile(pwd,'cropped_baseline_ambtemp');
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% ---------- Load meta CSV (preserve original headers) ----------
META = readtable(meta_csv, 'VariableNamingRule','preserve');

MouseNum = META.("Mouse");                 % numeric or string IDs: 1,2,10,...
Genotype = string(META.("Genotype"));
t_start  = META.("Time started (s)");
t_end    = META.("Time finished (s)");

OUT.entries = struct([]);
kEntry = 0;

for i = 1:numel(MouseNum)
    % mouse label used in filenames: "mouse2", "mouse10", etc.
    mouse_num   = MouseNum(i);
    mouse_label = "mouse" + string(mouse_num);
    geno        = Genotype(i);
    t0          = t_start(i);
    t1          = t_end(i);

    fprintf('\n=== Cropping %s (%s) from %.1f to %.1f s ===\n', ...
        mouse_label, geno, t0, t1);

    if isnan(t0) || isnan(t1) || t1 <= t0
        warning('  -> Invalid start/finish times, skipping this row.');
        continue;
    end

    % For each of the two conditions: baseline and ambtemp
    for cond = ["baseline","ambtemp"]

        fprintf('  Condition: %s\n', cond);

        % ---------- 1) Find scores CSV ----------
        % pattern: YYYYMMDD_<cond>_mouseX_GENO_scored_scores_1Hz.csv
        scores_pattern = sprintf('*_%s_%s_%s_scored_scores_1Hz.csv', ...
                                 cond, mouse_label, geno);
        scores_file = find_one_file(scores_dir, scores_pattern);
        if scores_file == ""
            warning('    No scores CSV found matching "%s".', scores_pattern);
            continue;
        end

        % ---------- 2) Find MAT file ----------
        % pattern: YYYYMMDD_<cond>_mouseX_GENO.mat
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

        keep = (T.time_s >= t0) & (T.time_s <= t1);
        T_crop = T(keep, :);
        % rebase time to 0 at manipulation start
        T_crop.time_s = T_crop.time_s - t0;

        [~, scores_base, ~] = fileparts(scores_file);
        out_scores_file = fullfile(out_dir, scores_base + "_crop.csv");
        writetable(T_crop, out_scores_file);
        fprintf('    -> Cropped scores: %s (n=%d rows)\n', out_scores_file, height(T_crop));

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
            i1 = min(nSamp, ceil(t1 * fs_eeg));

            eeg_crop = eeg(i0:i1);

            S_crop = S;         % start with original struct
            S_crop.(eegVar) = eeg_crop;
            S_crop.crop_t0_s = t0;
            S_crop.crop_t1_s = t1;
            S_crop.crop_cond = cond;
            S_crop.crop_mouse = mouse_label;
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
                j1 = min(nAmb, ceil(t1 * fs_amb));
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
        OUT.entries(kEntry).t_end_s     = t1;
        OUT.entries(kEntry).scores_in   = scores_file;
        OUT.entries(kEntry).scores_out  = out_scores_file;
        OUT.entries(kEntry).mat_in      = mat_file;
        OUT.entries(kEntry).mat_out     = out_mat_file;
    end
end

fprintf('\nDone. Cropped baseline + ambtemp segments are in: %s\n', out_dir);
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
