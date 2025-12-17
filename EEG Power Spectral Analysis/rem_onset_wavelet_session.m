function OUT = rem_onset_wavelet_session(eegFile, scoresCsv, outDir, varargin)
% rem_onset_wavelet_session
%
% For one baseline session:
%   - finds NREM→REM transitions from 1 Hz hypnogram
%   - extracts EEG around each REM onset (t = -t_pre .. +t_post)
%   - computes Morlet CWT for each segment (1–80 Hz)
%   - averages power across transitions -> time-freq map per mouse
%
% OUTPUT struct OUT (also saved to <outDir>/<mouse>_REM_onset_wavelet.mat):
%   OUT.F        : frequency vector (Hz)
%   OUT.tRel     : time relative to REM onset (s)
%   OUT.powMean  : mean power (freq x time, linear units)
%   OUT.pow_dB   : same, in dB (optional baseline norm)
%   OUT.meta     : mouseID, genotype, condition, nTransitions, etc.

    p = inputParser;
    addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
    addParameter(p,'epoch_sec',1,@(x)isscalar(x)&&x>0);
    addParameter(p,'t_pre',20,@(x)isscalar(x)&&x>0);   % seconds before REM onset
    addParameter(p,'t_post',40,@(x)isscalar(x)&&x>0);  % seconds after REM onset
    addParameter(p,'min_NREM_pre',10,@(x)isscalar(x)&&x>=0); % need this much pure NREM before
    addParameter(p,'min_REM_post',20,@(x)isscalar(x)&&x>=0); % and this much REM after
    parse(p,varargin{:});
    opt = p.Results;

    if nargin < 3 || isempty(outDir)
        outDir = fileparts(eegFile);
    end
    if ~exist(outDir,'dir'); mkdir(outDir); end

    % ---------- Load EEG ----------
    S = load(eegFile);
    if ~isfield(S,'eeg') || ~isfield(S,'eeg_frequency')
        error('Expected variables eeg and eeg_frequency in %s', eegFile);
    end
    eeg = double(S.eeg(:));
    fs  = double(S.eeg_frequency);

    % ---------- Load scores (1 Hz) ----------
    T = readtable(scoresCsv);

    % Adapt these lines to your CSV column names:
    % here assume a column "state" with numeric codes:
    if ismember('State', T.Properties.VariableNames)
        stateVec = T.state;
    else
        error('Adapt rem_onset_wavelet_session: could not find "state" column in %s', scoresCsv);
    end
    stateVec = double(stateVec(:));

    nEpoch = numel(stateVec);
    recDur = nEpoch * opt.epoch_sec;

    % Check EEG length
    t_eeg  = (0:numel(eeg)-1)/fs;
    if t_eeg(end) < recDur - opt.epoch_sec
        warning('EEG is shorter than score duration, clipping to EEG length.');
        nEpoch = floor(numel(eeg)/(fs*opt.epoch_sec));
        stateVec = stateVec(1:nEpoch);
        recDur   = nEpoch*opt.epoch_sec;
    end

    % ---------- Find NREM→REM transitions ----------
    c = opt.codes;
    isNREM = (stateVec == c.NREM);
    isREM  = (stateVec == c.REM);

    transIdx = find(isNREM(1:end-1) & isREM(2:end)); % index k: epoch k=NREM, k+1=REM
    fprintf('Found %d raw NREM→REM transitions.\n', numel(transIdx));

    % ---------- Filter transitions: need clean NREM before and REM after ----------
    t_pre_ep   = round(opt.t_pre / opt.epoch_sec);
    t_post_ep  = round(opt.t_post / opt.epoch_sec);
    minN_ep    = round(opt.min_NREM_pre / opt.epoch_sec);
    minR_ep    = round(opt.min_REM_post / opt.epoch_sec);

    keep = false(size(transIdx));
    for i = 1:numel(transIdx)
        k = transIdx(i);

        % indices in 1 Hz epochs relative to REM onset at epoch k+1
        kREM_on = k+1;

        % need window fully inside recording
        if (kREM_on - t_pre_ep) < 1 || (kREM_on + t_post_ep - 1) > nEpoch
            continue;
        end

        % Check NREM before
        preRange  = (kREM_on - minN_ep) : (kREM_on - 1);
        postRange = kREM_on : (kREM_on + minR_ep - 1);

        if all(stateVec(preRange) == c.NREM) && all(stateVec(postRange) == c.REM)
            keep(i) = true;
        end
    end
    transIdx = transIdx(keep);
    fprintf('Kept %d transitions after NREM/REM stability criteria.\n', numel(transIdx));

    if isempty(transIdx)
        warning('No usable NREM→REM transitions in this session.');
        OUT = [];
        return;
    end

    % ---------- Extract EEG segments and run wavelets ----------
    nTrans = numel(transIdx);
    t_pre_s  = opt.t_pre;
    t_post_s = opt.t_post;
    nSampWin = round((t_pre_s + t_post_s)*fs);

    % Time axis relative to REM onset
    tRel = linspace(-t_pre_s, t_post_s, nSampWin).';

    allPow = [];   % [nFreq x nTime x nTrans]

    for i = 1:nTrans
        k = transIdx(i);
        kREM_on = k+1;       % REM starts at epoch k+1

        % convert epoch index to sample index
        t0_s = (kREM_on*opt.epoch_sec) - t_pre_s;  % segment start (s)
        s0   = round(t0_s * fs) + 1;
        s1   = s0 + nSampWin - 1;

        if s0 < 1 || s1 > numel(eeg)
            continue; % safety check
        end

        seg = eeg(s0:s1);

        % Remove mean to avoid DC
        seg = seg - mean(seg);

        % Morlet-type CWT (analytic Morlet "amor"), 1–80 Hz
        [wt, f] = cwt(seg, fs, 'amor', 'FrequencyLimits', [1 80]);

        pow = abs(wt).^2;   % nFreq x nTime

        if isempty(allPow)
            allPow = zeros(size(pow,1), size(pow,2), nTrans);
        end
        allPow(:,:,i) = pow;
    end

    % ---------- Average across transitions ----------
    powMean = mean(allPow, 3, 'omitnan');

    % Optional: baseline-correct in dB relative to pre-REM NREM
    baseMask = tRel < 0;  % all pre-onset
    basePow  = mean(powMean(:, baseMask), 2, 'omitnan');  % freq x 1
    pow_dB   = 10*log10(bsxfun(@rdivide, powMean, basePow)); % relative change in dB

    % ---------- Meta (try to parse from filename) ----------
    [~, baseName] = fileparts(eegFile);
    meta = struct();
    meta.eegFile   = eegFile;
    meta.scoresCsv = scoresCsv;
    meta.mouseID   = baseName;  % you can override with your own parser
    meta.nTransitions = nTrans;
    meta.fs        = fs;
    meta.t_pre     = t_pre_s;
    meta.t_post    = t_post_s;

    OUT = struct();
    OUT.F       = f;
    OUT.tRel    = tRel;
    OUT.powMean = powMean;
    OUT.pow_dB  = pow_dB;
    OUT.meta    = meta;

    % ---------- Save ----------
    [~, nameOnly] = fileparts(eegFile);
    outFile = fullfile(outDir, sprintf('%s_REM_onset_wavelet.mat', nameOnly));
    save(outFile, 'OUT');
    fprintf('Saved REM-onset wavelet to %s\n', outFile);
end
