function select_REM_for_TFR_context(varargin)
% Build a REM-selection .mat ready for rem_timefreq_analysis:
%  - bouts(k).EEG (and EMG if present) contain CONTEXT-padded segments
%  - rem_bouts_table reports CORE REM (no padding) for correct bout metrics
%  - top-level eeg_frequency (and eeg) included for convenience
%
% Usage:
% select_REM_for_TFR_context('file','your_scored.mat','PRE_SEC',3,'POST_SEC',3);

%% -------- Args --------
p = inputParser;
addParameter(p,'file','',@(s)ischar(s)||isstring(s));
addParameter(p,'VARS',struct('eeg','eeg','emg','emg','scores','sleep_scores', ...
                             'fs_eeg','eeg_frequency'));
addParameter(p,'EPOCH_SEC',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'REM_CODES',[2],@(v)isnumeric(v)&&~isempty(v));
addParameter(p,'OUTDIR','',@(s)ischar(s)||isstring(s));
addParameter(p,'PRE_SEC',3,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'POST_SEC',3,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
parse(p,varargin{:});

FILE=string(p.Results.file); VARS=p.Results.VARS; EPOCH_SEC=p.Results.EPOCH_SEC;
REM_CODES=p.Results.REM_CODES; OUTDIR=string(p.Results.OUTDIR);
PRE_SEC=p.Results.PRE_SEC; POST_SEC=p.Results.POST_SEC;

if FILE=="" || ~isfile(FILE)
    [f,fp]=uigetfile('*.mat','Select scored .mat'); assert(f~=0); FILE=string(fullfile(fp,f));
end
if OUTDIR=="", OUTDIR=string(fileparts(FILE)); end
%% -------- Load --------
S = load(FILE);
assert(isfield(S,VARS.scores),'Missing scores "%s".',VARS.scores);
assert(isfield(S,VARS.eeg),'Missing EEG "%s".',VARS.eeg);
scores = S.(VARS.scores)(:);
eeg    = S.(VARS.eeg)(:);
hasEMG = isfield(S,VARS.emg);
if hasEMG, emg = S.(VARS.emg)(:); else, emg = []; end

fs_eeg = []; if isfield(S,VARS.fs_eeg), fs_eeg = double(S.(VARS.fs_eeg)); end
assert(~isempty(fs_eeg)&&isscalar(fs_eeg),'Missing/invalid "%s".',VARS.fs_eeg);

%% -------- Core REM bouts (metrics) --------
nS = numel(eeg);
t  = (0:nS-1)'/fs_eeg;
samp_per_ep = round(EPOCH_SEC*fs_eeg);

isREM_ep = ismember(scores, REM_CODES);
d = diff([false; isREM_ep; false]);
start_ep = find(d==1); end_ep = find(d==-1)-1;
nB = numel(start_ep);

core_s0 = (start_ep-1)*samp_per_ep + 1;
core_s1 = min(end_ep*samp_per_ep, nS);

% ------- table for metrics (core only) -------
start_time_s = (start_ep - 1) * EPOCH_SEC;
end_time_s   =  end_ep * EPOCH_SEC;
dur_s        = (end_ep - start_ep + 1) * EPOCH_SEC;

rem_bouts_table = table( ...
    start_ep(:), end_ep(:), ...
    start_time_s(:), end_time_s(:), dur_s(:), ...
    'VariableNames', {'start_epoch','end_epoch','start_time_s','end_time_s','dur_s'});

% ------- per-sample core mask (for completeness) -------
rem_mask_samples = false(nS,1);
for k = 1:nB
    s0 = core_s0(k); 
    s1 = core_s1(k);
    if s0 <= s1
        rem_mask_samples(s0:s1) = true;
    end
end


%% -------- Build CONTEXT bouts (for TFR input) --------
preN  = round(PRE_SEC*fs_eeg);
postN = round(POST_SEC*fs_eeg);

bouts = struct('idx_samples',[],'t_s',[],'EEG',[],'EMG',[]);
bouts = repmat(bouts, nB, 1);

for k=1:nB
    s0p = max(1, core_s0(k)-preN);
    s1p = min(nS, core_s1(k)+postN);
    idx = (s0p:s1p)';
    bouts(k).idx_samples = idx;
    bouts(k).t_s         = t(idx);
    bouts(k).EEG         = eeg(idx);
    if hasEMG, bouts(k).EMG = emg(idx); end
end

%% -------- Params & Save --------
params = struct('file',FILE,'EPOCH_SEC',EPOCH_SEC,'REM_CODES',REM_CODES, ...
                'fs_eeg',fs_eeg,'PRE_SEC',PRE_SEC,'POST_SEC',POST_SEC);

[~,base] = fileparts(FILE);
out_mat  = fullfile(OUTDIR, base + "_REMselection_CONTEXT_forTFR.mat");

% Fields rem_timefreq_analysis looks for:
% - eeg_frequency (fs), bouts, rem_bouts_table; having 'eeg' at top is a plus for 'zrem'
eeg_frequency = fs_eeg; 
save(out_mat, 'eeg','eeg_frequency','bouts','rem_bouts_table','rem_mask_samples','params','-v7.3');

fprintf('[OK] Wrote %s\n', out_mat);
fprintf('  bouts: %d  (core metrics unchanged; segments include %.1fs pre / %.1fs post)\n', nB, PRE_SEC, POST_SEC);
end
