function OUT = run_group_psd_APPvsWT(eegDir, scoreDir, out_dir, varargin)
% run_group_psd_APPvsWT
% -------------------------------------------------------------------------
% Compute group Welch PSDs for WT vs APP, per state, for ONE condition
% (e.g. baseline or ambtemp) and plot:
%
%   3 subplots: WAKE, NREM, REM
%   - grey line + shaded SEM  : WT
%   - blue  line + shaded SEM : APP
%
% ASSUMPTIONS
%   - EEG .mat filenames: YYYYMMDD_condition_mouseX_GENO.mat
%       e.g. 20251003_ambtemp_mouse1_APP.mat
%   - Matching scores CSV: YYYYMMDD_condition_mouseX_GENO_scored_scores_1Hz.csv
%   - .mat has variables: eeg, eeg_frequency
%   - read_scores_csv(scores_csv) returns:
%       [states, epochs_t, ~, mouse_id, geno, cond]
%   - states codes (change if needed):
%       WK=0, NREM=1, REM=2
%
% NAME–VALUE OPTIONS
%   'condition'   : which condition to include (default 'ambtemp')
%   'epoch_sec'   : length of each scoring epoch (default 1 s)
%   'minREMsec'   : minimum total seconds per state (default 60)
%   'nfft'        : NFFT for pwelch (default 1024)
%   'welch_win_s' : window length in seconds (default 4)
%   'welch_olap'  : fraction overlap, 0–1 (default 0.5)
%
% OUTPUT
%   OUT.f_Hz                 : frequency vector
%   OUT.psd_WT.(State)       : [nWT x nFreq] PSD (dB)
%   OUT.psd_APP.(State)      : [nAPP x nFreq] PSD (dB)
%   OUT.fig_file             : path to the PNG figure
% -------------------------------------------------------------------------

p = inputParser;
addParameter(p,'condition','ambtemp',@(s)ischar(s)||isstring(s));
addParameter(p,'epoch_sec',1,@(x)isscalar(x)&&x>0);
addParameter(p,'minREMsec',60,@(x)isscalar(x)&&x>=0);
addParameter(p,'nfft',1024,@(x)isscalar(x)&&x>0);
addParameter(p,'welch_win_s',4,@(x)isscalar(x)&&x>0);
addParameter(p,'welch_olap',0.5,@(x)isscalar(x)&&x>=0&&x<1);
parse(p,varargin{:});

condWanted  = lower(strtrim(string(p.Results.condition)));
epoch_sec   = p.Results.epoch_sec;
minStateSec = p.Results.minREMsec;
NFFT        = p.Results.nfft;
win_s       = p.Results.welch_win_s;
olap_frac   = p.Results.welch_olap;

if nargin < 3 || isempty(out_dir)
    out_dir = fullfile(pwd, sprintf('EEG_PSD_%s',condWanted));
end
if ~isfolder(out_dir), mkdir(out_dir); end

% ---------- State code mapping (EDIT if your codes differ) ----------
codes.WAKE = 0;
codes.NREM = 1;
codes.REM  = 2;

statesWanted = {'WAKE','NREM','REM'};

COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

files = dir(fullfile(eegDir, '*.mat'));
assert(~isempty(files), 'No .mat files found in %s', eegDir);

% Containers for PSDs
psd_WT  = struct();   % fields WAKE / NREM / REM -> [nWT x nFreq]
psd_APP = struct();
for s = 1:numel(statesWanted)
    st = statesWanted{s};
    psd_WT.(st)  = [];
    psd_APP.(st) = [];
end

f_common = [];
fprintf('=== Group PSD for condition: %s ===\n', condWanted);

for k = 1:numel(files)
    mat_path = fullfile(files(k).folder, files(k).name);
    [~, base, ~] = fileparts(files(k).name);

    % Parse filename: YYYYMMDD_cond_mouseX_GENO
    tok = regexp(base, ...
        '^(?<date>\d{8})_(?<cond>[^_]+)_(?<mouse>mouse\d+)_(?<geno>APP|WT)$', ...
        'names');

    if isempty(tok)
        fprintf('  Skipping %s (name not matching pattern)\n', base);
        continue;
    end

    cond = lower(string(tok.cond));
    if cond ~= condWanted
        continue;   % wrong condition
    end

    mouseID = string(tok.mouse);
    geno    = string(tok.geno);

    % Matching scores CSV name
    scores_name = sprintf('%s_%s_%s_%s_scored_scores_1Hz.csv', ...
                          tok.date, tok.cond, tok.mouse, tok.geno);
    scores_path = fullfile(scoreDir, scores_name);
    if ~isfile(scores_path)
        fprintf('  No scores CSV for %s, expected %s, skipping.\n', ...
            base, scores_name);
        continue;
    end

    fprintf('  %s (%s, %s) -> %s\n', base, geno, cond, scores_name);

    % ---- Load EEG ----
    S = load(mat_path);
    if isfield(S,'eeg_notched')
        eeg = S.eeg_notched;
    elseif isfield(S,'eeg')
        eeg = S.eeg;
    else
        warning('File %s has no eeg/eeg_notched variable, skipping.', base);
        continue;
    end
    if isfield(S,'eeg_frequency')
        fs = S.eeg_frequency;
    else
        error('File %s missing eeg_frequency.', base);
    end
    eeg = double(eeg(:));

    % ---- Scores ----
    [states, epochs_t, ~, ~, ~, ~] = read_scores_csv(scores_path);
    states = double(states(:));
    nEpoch = numel(states);

    % Ensure epochs mapping to EEG length
    samples_per_epoch = round(epoch_sec * fs);
    if numel(eeg) < nEpoch * samples_per_epoch
        warning('EEG shorter than states*epoch_sec for %s, truncating states.', base);
        nEpoch = floor(numel(eeg)/samples_per_epoch);
        states = states(1:nEpoch);
    end

    % reshape EEG into [samplesPerEpoch x nEpoch]
    eeg = eeg(1:nEpoch*samples_per_epoch);
    eeg_mat = reshape(eeg, samples_per_epoch, nEpoch);  % samples x epochs

    % Welch parameters
    win  = hamming(round(win_s*fs));
    nover= round(olap_frac * numel(win));

    % ---- Per state ----
    for s = 1:numel(statesWanted)
        st = statesWanted{s};
        code = codes.(st);

        mask_epochs = (states == code);
        nE = sum(mask_epochs);
        if nE==0, continue; end

        % concatenate all epochs of this state
        seg = eeg_mat(:, mask_epochs);
        seg = seg(:);
        totalSec = numel(seg)/fs;

        if totalSec < minStateSec
            % too little data -> noisy, skip
            continue;
        end

        % Welch PSD
        [Pxx, f] = pwelch(seg, win, nover, NFFT, fs);
        Pxx_dB = 10*log10(Pxx);

        if isempty(f_common)
            f_common = f(:);
        else
            % ensure same freqs; if not, interp
            if numel(f) ~= numel(f_common) || any(abs(f-f_common)>1e-6)
                Pxx_dB = interp1(f, Pxx_dB, f_common, 'linear', 'extrap');
            end
        end

        % append to genotype-specific matrix
        if geno == "WT"
            psd_WT.(st) = [psd_WT.(st); Pxx_dB(:).'];
        else
            psd_APP.(st) = [psd_APP.(st); Pxx_dB(:).'];
        end
    end
end

if isempty(f_common)
    warning('No PSDs computed – check condition / file patterns.');
    OUT = struct('success',false);
    return;
end

% ---------- Plot group figure ----------
fig = figure('Color','w','Position',[100 100 1400 450]);

for s = 1:numel(statesWanted)
    st = statesWanted{s};
    subplot(1, numel(statesWanted), s); hold on;

    WTmat  = psd_WT.(st);
    APPmat = psd_APP.(st);

    haveWT  = ~isempty(WTmat);
    haveAPP = ~isempty(APPmat);

    if ~haveWT && ~haveAPP
        title(st);
        continue;
    end

    if haveWT
        mWT  = mean(WTmat,1,'omitnan');
        seWT = std(WTmat,[],1,'omitnan') ./ sqrt(size(WTmat,1));
        % shaded area
        fill([f_common; flipud(f_common)], ...
             [mWT-seWT, fliplr(mWT+seWT)], ...
             COL_WT, 'FaceAlpha',0.3,'EdgeColor','none');
        plot(f_common, mWT, 'Color',COL_WT, 'LineWidth',1.5);
    end

    if haveAPP
        mAPP  = mean(APPmat,1,'omitnan');
        seAPP = std(APPmat,[],1,'omitnan') ./ sqrt(size(APPmat,1));
        fill([f_common; flipud(f_common)], ...
             [mAPP-seAPP, fliplr(mAPP+seAPP)], ...
             COL_APP, 'FaceAlpha',0.3,'EdgeColor','none');
        plot(f_common, mAPP, 'Color',COL_APP, 'LineWidth',1.5);
    end

    xlim([0 100]);
    ylim([min([WTmat(:); APPmat(:)])-5, max([WTmat(:); APPmat(:)])+5]);
    ylabel('Power (dB)');
    xlabel('Frequency (Hz)');
    title(st, 'Interpreter','none');
    set(gca,'Box','off');

    if s==1
        legend({'WT SEM','WT mean','APP SEM','APP mean'}, ...
               'Location','southwest');
    end
end

sgtitle(sprintf('Condition: %s', condWanted), 'FontWeight','bold');

fig_file = fullfile(out_dir, sprintf('PSD_%s_APPvsWT.png', condWanted));
saveas(fig, fig_file);

OUT = struct();
OUT.success   = true;
OUT.f_Hz      = f_common;
OUT.psd_WT    = psd_WT;
OUT.psd_APP   = psd_APP;
OUT.fig_file  = fig_file;

fprintf('✅ Group PSD figure saved to %s\n', fig_file);
end
