function compare_emg_psd_ambtemp(sigDir, maxDurMin)
% COMPARE_EMG_PSD_AMBTEMP
% -------------------------------------------------------------------------
% Scan a folder for EMG .mat files, select only "ambtemp" recordings, and
% plot their EMG power spectral densities (PSDs) to help identify
% artifact frequencies (e.g. warming lamp).
%
% Assumes filenames like:
%   YYYYMMDD-cond-mouseNN_notched50Hz_bw3.mat
% and selects only those with cond = 'ambtemp'.
%
% INPUTS (optional)
%   sigDir    : folder with .mat signal files
%               (if empty or omitted -> dialog to choose)
%   maxDurMin : max duration (in minutes) from each recording used
%               for PSD (default = 10 min)
%
% OUTPUT
%   Just plots a figure with all PSDs overlaid + their mean.
% -------------------------------------------------------------------------

if nargin < 1 || isempty(sigDir)
    sigDir = uigetdir(pwd, 'Select folder with EMG .mat signal files');
    if sigDir == 0
        error('No folder selected.');
    end
end

if nargin < 2 || isempty(maxDurMin)
    maxDurMin = 10;  % use first 10 minutes by default
end

fprintf('\n[compare_emg_psd_ambtemp]\nSignal dir: %s\n', sigDir);

% --------------------- 1) find ambtemp .mat files -----------------------
allMat = dir(fullfile(sigDir, '*.mat'));
keep   = false(numel(allMat),1);
meta   = struct([]);

for k = 1:numel(allMat)
    fn = allMat(k).name;
    info = parse_signal_filename_cond(fn);
    if isempty(info), continue; end

    if strcmpi(info.cond, 'ambtemp')
        keep(k) = true;

        info.file = fullfile(sigDir, fn);
        if isempty(meta)
            meta = info;
        else
            meta(end+1) = info; %#ok<AGROW>
        end
    end
end

if isempty(meta)
    error('No "ambtemp" .mat files found in folder %s', sigDir);
end

allMat = allMat(keep);
nFiles = numel(meta);
fprintf('Found %d ambtemp recordings.\n', nFiles);

% --------------------- 2) loop over ambtemp files, compute PSD ----------
PSD_all = {};
F_all   = {};
labels  = cell(nFiles,1);

for k = 1:nFiles
    fn = meta(k).file;
    fprintf('\n[%d/%d] %s\n', k, nFiles, fn);

    % ---- load EMG + fs (robust to different variable names) ----
    S = load(fn);

    cand_emg = {'emg','EMG'};
    cand_fs  = {'emg_frequency','eeg_frequency','Fs_emg','fs_emg'};

    emg = [];
    for i = 1:numel(cand_emg)
        if isfield(S, cand_emg{i})
            emg = S.(cand_emg{i});
            emg_name = cand_emg{i};
            break;
        end
    end
    if isempty(emg)
        fprintf('  [WARN] No EMG variable found in %s, skipping.\n', fn);
        continue;
    end
    emg = double(emg(:));

    fs_emg = [];
    for i = 1:numel(cand_fs)
        if isfield(S, cand_fs{i})
            fs_emg = S.(cand_fs{i});
            fs_name = cand_fs{i};
            break;
        end
    end
    if isempty(fs_emg)
        fprintf('  [WARN] No EMG sampling rate found in %s, skipping.\n', fn);
        continue;
    end
    fs_emg = double(fs_emg);

    fprintf('  Using EMG: %s, fs = %.2f Hz (%s)\n', emg_name, fs_emg, fs_name);

    % ---- take first maxDurMin minutes (or full length if shorter) ----
    maxSamples = round(maxDurMin * 60 * fs_emg);
    if numel(emg) > maxSamples
        emg_seg = emg(1:maxSamples);
    else
        emg_seg = emg;
    end

    % ---- compute PSD using Welch ----
    % 4-second window, 50% overlap, full frequency range up to fs/2
    win   = round(4 * fs_emg);
    nover = round(win / 2);
    nfft  = max(2^nextpow2(win), 512);

    [Pxx, F] = pwelch(emg_seg, win, nover, nfft, fs_emg);

    PSD_all{end+1} = Pxx; %#ok<AGROW>
    F_all{end+1}   = F;
    labels{k}      = sprintf('mouse%s (%s)', meta(k).mouse, meta(k).date);
end

% Remove any empties (files that were skipped)
valid_idx = ~cellfun(@isempty, PSD_all);
PSD_all = PSD_all(valid_idx);
F_all   = F_all(valid_idx);
labels  = labels(valid_idx);

if isempty(PSD_all)
    error('No valid EMG PSDs computed (missing EMG or fs in files).');
end

% --------------------- 3) unify frequency vector if needed --------------
% In practice, if fs_emg is same for all, F vectors will match.
F_ref = F_all{1};
nF    = numel(F_ref);

PSD_mat = nan(nF, numel(PSD_all));
for k = 1:numel(PSD_all)
    if numel(F_all{k}) == nF && all(abs(F_all{k} - F_ref) < 1e-6)
        PSD_mat(:,k) = PSD_all{k};
    else
        % interpolate onto reference grid
        PSD_mat(:,k) = interp1(F_all{k}, PSD_all{k}, F_ref, 'linear', 'extrap');
    end
end

% Convert to dB for plotting
PSD_dB = 10*log10(PSD_mat);

% --------------------- 4) plot all PSDs + mean --------------------------
figure('Name','EMG PSD - ambtemp recordings','Color','w');

subplot(1,2,1); hold on;
for k = 1:size(PSD_dB,2)
    plot(F_ref, PSD_dB(:,k), 'Color',[0.7 0.7 0.7]);
end
plot(F_ref, mean(PSD_dB,2,'omitnan'), 'k','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('EMG PSD (first %.1f min) - all ambtemp', maxDurMin));
xlim([0 max(F_ref)]);
box off;
legend([repmat({''},1,size(PSD_dB,2)) {'mean'}], 'Location','best'); % simple legend

% Zoom into 0–200 Hz with labels
subplot(1,2,2); hold on;
for k = 1:size(PSD_dB,2)
    plot(F_ref, PSD_dB(:,k), 'Color',[0.7 0.7 0.7]);
end
plot(F_ref, mean(PSD_dB,2,'omitnan'), 'k','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Zoom 0–200 Hz');
xlim([0 200]);
box off;

% Annotate with number of recordings
sgtitle(sprintf('EMG PSD in ambtemp (%d recordings)', size(PSD_dB,2)));

% Also print out a little summary to console to help you read peaks
fprintf('\nApproximate peak frequencies per recording (top ~3 peaks under 200 Hz):\n');
for k = 1:size(PSD_dB,2)
    [pks, locs] = findpeaks(PSD_dB(:,k), F_ref, 'MinPeakProminence', 5);
    sel = locs(locs <= 200);
    psel = pks(locs <= 200);
    [~, idxSort] = sort(psel, 'descend');
    idxSort = idxSort(1:min(3,numel(idxSort)));
    fprintf('  %s: ', labels{k});
    fprintf('%.1f Hz  ', sel(idxSort));
    fprintf('\n');
end

end

% =======================================================================
function info = parse_signal_filename_cond(fname)
% Parse: YYYYMMDD-cond-mouseNN_notched50Hz_bw3.mat
% Returns empty [] if it doesn't match.
pat = '^(?<date>\d{8})-(?<cond>[^-]+)-mouse(?<mouse>\d+).*\.mat$';
m = regexp(fname, pat, 'names');
if isempty(m)
    info = [];
else
    info = struct();
    info.date  = m.date;
    info.cond  = lower(m.cond);
    info.mouse = m.mouse;
end
end
