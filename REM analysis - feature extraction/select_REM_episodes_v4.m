function select_REM_episodes_v4(varargin)
% REM selection with strict vs padded exports to separate folders
% Preserves core bout metrics; padding only affects extracted segments.

% -------------------- Args --------------------
p = inputParser;
addParameter(p,'file','',@(s)ischar(s)||isstring(s));
addParameter(p,'VARS',struct('eeg','eeg','emg','emg','ach','ne', ...
                             'scores','sleep_scores','fs_eeg','eeg_frequency','fs_ach','ne_frequency'));
addParameter(p,'EPOCH_SEC',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'REM_CODES',[2],@(v)isnumeric(v)&&~isempty(v));
addParameter(p,'OUTDIR','',@(s)ischar(s)||isstring(s));

% padding
addParameter(p,'PRE_SEC',3,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'POST_SEC',3,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MERGE_PADDED',false,@islogical);

% NEW: separate export toggles/paths
addParameter(p,'EXPORT_SEPARATE',true,@islogical);
addParameter(p,'STRICT_DIR','REM_STRICT',@(s)ischar(s)||isstring(s));
addParameter(p,'PAD_DIR','REM_PADDED',@(s)ischar(s)||isstring(s));

parse(p,varargin{:});
FILE       = string(p.Results.file);
VARS       = p.Results.VARS;
EPOCH_SEC  = p.Results.EPOCH_SEC;
REM_CODES  = p.Results.REM_CODES;
OUTDIR     = string(p.Results.OUTDIR);
PRE_SEC    = p.Results.PRE_SEC;
POST_SEC   = p.Results.POST_SEC;
MERGE_PAD  = p.Results.MERGE_PADDED;
EXPORT_SEP = p.Results.EXPORT_SEPARATE;
STRICT_DIR = string(p.Results.STRICT_DIR);
PAD_DIR    = string(p.Results.PAD_DIR);

if FILE=="" || ~isfile(FILE)
    [f,fp] = uigetfile('*.mat','Select scored .mat');
    assert(f~=0,'No file selected.');
    FILE = string(fullfile(fp,f));
end
if OUTDIR=="", OUTDIR = string(fileparts(FILE)); end

% -------------------- Load --------------------
S = load(FILE);
assert(isfield(S,VARS.scores), 'Missing scores variable "%s".', VARS.scores);
scores = S.(VARS.scores)(:);

assert(isfield(S,VARS.eeg), 'Missing EEG variable "%s".', VARS.eeg);
EEG = S.(VARS.eeg)(:);

fs_eeg = []; fs_ach = [];
if isfield(S,VARS.fs_eeg), fs_eeg = S.(VARS.fs_eeg); end
assert(~isempty(fs_eeg) && isnumeric(fs_eeg) && isscalar(fs_eeg), ...
    'Missing/invalid fs_eeg "%s".', VARS.fs_eeg);

hasEMG = isfield(S,VARS.emg);   if hasEMG, EMG = S.(VARS.emg)(:); else, EMG = []; end
hasACH = isfield(S,VARS.ach);   if hasACH, ACH = S.(VARS.ach)(:); else, ACH = []; end
if hasACH && isfield(S,VARS.fs_ach), fs_ach = S.(VARS.fs_ach); end

% -------------------- Epoch sanity --------------------
nSamples       = numel(EEG);
total_dur_s    = nSamples/fs_eeg;
expected_ep    = floor(total_dur_s/EPOCH_SEC);
if abs(numel(scores)-expected_ep) > max(5,0.01*expected_ep)
    warning(['scores length (%d) differs from expected (%d). ' ...
             'Check EPOCH_SEC=%.3f and fs_eeg=%.3f Hz.'], numel(scores), expected_ep, EPOCH_SEC, fs_eeg);
end

% -------------------- Core REM bouts (epoch space) --------------------
isREM_ep = ismember(scores, REM_CODES);   % core only
d = diff([false; isREM_ep; false]);
start_epoch = find(d==1);
end_epoch   = find(d==-1)-1;
nBouts      = numel(start_epoch);

samples_per_epoch = round(EPOCH_SEC*fs_eeg);
t_eeg = (0:nSamples-1)'/fs_eeg;

% Core boundaries
core_start_sample = (start_epoch-1)*samples_per_epoch + 1;
core_end_sample   = min(end_epoch*samples_per_epoch, nSamples);
core_start_time_s = (start_epoch-1)*EPOCH_SEC;
core_end_time_s   = end_epoch*EPOCH_SEC;
core_dur_s        = (end_epoch - start_epoch + 1)*EPOCH_SEC;

rem_bouts_table_core = table( ...
    start_epoch(:), end_epoch(:), ...
    core_start_sample(:), core_end_sample(:), ...
    core_start_time_s(:), core_end_time_s(:), core_dur_s(:), ...
    'VariableNames',{'start_epoch','end_epoch','core_start_sample','core_end_sample', ...
                     'core_start_time_s','core_end_time_s','core_dur_s'});

% Core masks/IDs (for metrics)
rem_mask_samples_core   = false(nSamples,1);
bout_id_per_sample_core = zeros(nSamples,1,'uint32');
for k=1:nBouts
    if core_start_sample(k)<=core_end_sample(k)
        rem_mask_samples_core(core_start_sample(k):core_end_sample(k)) = true;
        bout_id_per_sample_core(core_start_sample(k):core_end_sample(k)) = k;
    end
end

% -------------------- Padded boundaries (for extraction only) -----------
pre_samp  = round(PRE_SEC*fs_eeg);
post_samp = round(POST_SEC*fs_eeg);

pad_start_sample = max(1, core_start_sample - pre_samp);
pad_end_sample   = min(nSamples, core_end_sample + post_samp);

pad_start_time_s = (pad_start_sample-1)/fs_eeg;
pad_end_time_s   = (pad_end_sample-1)/fs_eeg;
pad_dur_s        = pad_end_time_s - pad_start_time_s;

rem_bouts_table_padded = table( ...
    pad_start_sample(:), pad_end_sample(:), ...
    pad_start_time_s(:), pad_end_time_s(:), pad_dur_s(:), ...
    'VariableNames',{'pad_start_sample','pad_end_sample','pad_start_time_s','pad_end_time_s','pad_dur_s'});

% Optional merged padded windows
if p.Results.MERGE_PADDED && nBouts>0
    A = [pad_start_sample(:), pad_end_sample(:)];
    A = sortrows(A,1);
    merged = [];
    cur = A(1,:);
    for i=2:size(A,1)
        if A(i,1) <= cur(2)+1
            cur(2) = max(cur(2), A(i,2));
        else
            merged = [merged; cur]; %#ok<AGROW>
            cur = A(i,:);
        end
    end
    pad_merged = [merged; cur];
else
    pad_merged = [];
end

% -------------------- Extract segments --------------------
bouts = struct('idx_core',[],'idx_padded',[],'t_core',[],'t_padded',[], ...
               'EEG_core',[],'EEG_padded',[],'EMG_core',[],'EMG_padded',[], ...
               'ACH_core',[],'ACH_padded',[]);
bouts = repmat(bouts, nBouts, 1);

for k=1:nBouts
    % Core
    ic = (core_start_sample(k):core_end_sample(k))';
    bouts(k).idx_core = ic;
    bouts(k).t_core   = t_eeg(ic);
    bouts(k).EEG_core = EEG(ic);
    if hasEMG, bouts(k).EMG_core = EMG(ic); end
    if hasACH && ~isempty(fs_ach)
        tc = (ic-1)/fs_eeg;
        ia = round(tc*fs_ach)+1; ia = min(max(ia,1), numel(ACH));
        bouts(k).ACH_core = ACH(ia);
    end

    % Padded
    ip = (pad_start_sample(k):pad_end_sample(k))';
    bouts(k).idx_padded = ip;
    bouts(k).t_padded   = t_eeg(ip);
    bouts(k).EEG_padded = EEG(ip);
    if hasEMG, bouts(k).EMG_padded = EMG(ip); end
    if hasACH && ~isempty(fs_ach)
        tp = (ip-1)/fs_eeg;
        ia = round(tp*fs_ach)+1; ia = min(max(ia,1), numel(ACH));
        bouts(k).ACH_padded = ACH(ia);
    end
end

% Optional merged padded segments
merged_segments = struct('idx',[],'t',[],'EEG',[],'EMG',[],'ACH',[]);
if ~isempty(pad_merged)
    merged_segments = repmat(merged_segments, size(pad_merged,1), 1);
    for i=1:size(pad_merged,1)
        ii = (pad_merged(i,1):pad_merged(i,2))';
        merged_segments(i).idx = ii;
        merged_segments(i).t   = t_eeg(ii);
        merged_segments(i).EEG = EEG(ii);
        if hasEMG, merged_segments(i).EMG = EMG(ii); end
        if hasACH && ~isempty(fs_ach)
            tt = (ii-1)/fs_eeg;
            ia = round(tt*fs_ach)+1; ia = min(max(ia,1), numel(ACH));
            merged_segments(i).ACH = ACH(ia);
        end
    end
end

% -------------------- Save master bundle (as before) --------------------
params = struct('file',FILE,'VARS',VARS,'EPOCH_SEC',EPOCH_SEC,'REM_CODES',REM_CODES, ...
                'fs_eeg',fs_eeg,'fs_ach',fs_ach, ...
                'PRE_SEC',PRE_SEC,'POST_SEC',POST_SEC,'MERGE_PADDED',MERGE_PAD, ...
                'nSamples',nSamples,'nEpochs',numel(scores),'nBouts',nBouts);

[~, base] = fileparts(FILE);
out_mat_master = fullfile(OUTDIR, base + "_REMselection_PAD.mat");
out_csv_core   = fullfile(OUTDIR, base + "_REMbouts_core.csv");
out_csv_padded = fullfile(OUTDIR, base + "_REMbouts_padded.csv");

save_vars = S; % keep original variables
save_vars.rem_mask_samples_core    = rem_mask_samples_core;
save_vars.bout_id_per_sample_core  = bout_id_per_sample_core;
save_vars.rem_bouts_table_core     = rem_bouts_table_core;
save_vars.rem_bouts_table_padded   = rem_bouts_table_padded;
if ~isempty(pad_merged)
    rem_padded_merged_table = table( ...
        pad_merged(:,1), pad_merged(:,2), ...
        (pad_merged(:,1)-1)/fs_eeg, (pad_merged(:,2)-1)/fs_eeg, ...
        ((pad_merged(:,2)-pad_merged(:,1))/fs_eeg), ...
        'VariableNames',{'pad_start_sample','pad_end_sample','pad_start_time_s','pad_end_time_s','pad_dur_s'});
    save_vars.rem_padded_merged_table = rem_padded_merged_table;
    save_vars.padded_merged_segments  = merged_segments;
end
save_vars.bouts_with_padding = bouts;
save_vars.params_REMselection_PAD = params;
save(out_mat_master,'-struct','save_vars','-v7.3');
try, writetable(rem_bouts_table_core, out_csv_core);   catch, end
try, writetable(rem_bouts_table_padded, out_csv_padded);catch, end

fprintf('[OK] Master saved:\n  %s\n  %s\n  %s\n', out_mat_master, out_csv_core, out_csv_padded);
fprintf('[NOTE] Bout counts/durations use CORE only.\n');

% -------------------- NEW: Separate strict/padded exports ---------------
if EXPORT_SEP
    strict_root = fullfile(OUTDIR, base + "__" + STRICT_DIR);
    padded_root = fullfile(OUTDIR, base + "__" + PAD_DIR);
    if ~exist(strict_root,'dir'), mkdir(strict_root); end
    if ~exist(padded_root,'dir'), mkdir(padded_root); end

    % STRICT combined + CSV
    strict_combined = struct('bouts_core',[],'fs_eeg',fs_eeg,'fs_ach',fs_ach,'VARS',VARS,'params',params);

    strict_combined.bouts_core = arrayfun(@(b) struct( ...
        'idx', b.idx_core, ...
        't',   b.t_core, ...
        'EEG', b.EEG_core, ...
        'EMG', b.EMG_core, ...   % <- no getfield; these exist (possibly [])
        'ACH', b.ACH_core), ...
        bouts, 'UniformOutput', false);

    save(fullfile(strict_root, 'REM_STRICT_combined.mat'),'-struct','strict_combined','-v7.3');
    writetable(rem_bouts_table_core, fullfile(strict_root,'REM_STRICT_bouts.csv'));

    % STRICT per-bout
    for k=1:nBouts
        B = struct('idx',bouts(k).idx_core,'t',bouts(k).t_core,'EEG',bouts(k).EEG_core);
        if hasEMG, B.EMG = bouts(k).EMG_core; end
        if hasACH, B.ACH = bouts(k).ACH_core; end
        save(fullfile(strict_root, sprintf('bout_%03d.mat',k)),'-struct','B','-v7.3');
    end

    % PADDED combined + CSV
    % PADDED combined + CSV
    padded_combined = struct('bouts_padded',[],'fs_eeg',fs_eeg,'fs_ach',fs_ach,'VARS',VARS,'params',params);

    padded_combined.bouts_padded = arrayfun(@(b) struct( ...
        'idx', b.idx_padded, ...
        't',   b.t_padded, ...
        'EEG', b.EEG_padded, ...
        'EMG', b.EMG_padded, ...
        'ACH', b.ACH_padded), ...
        bouts, 'UniformOutput', false);

    save(fullfile(padded_root, 'REM_PADDED_combined.mat'),'-struct','padded_combined','-v7.3');
    writetable(rem_bouts_table_padded, fullfile(padded_root,'REM_PADDED_bouts.csv'));


    % PADDED per-bout
    for k=1<nBouts
        B = struct('idx',bouts(k).idx_padded,'t',bouts(k).t_padded,'EEG',bouts(k).EEG_padded);
        if hasEMG, B.EMG = bouts(k).EMG_padded; end
        if hasACH, B.ACH = bouts(k).ACH_padded; end
        save(fullfile(padded_root, sprintf('bout_%03d.mat',k)),'-struct','B','-v7.3');
    end

    fprintf('[OK] Separate exports written to:\n  %s\n  %s\n', strict_root, padded_root);
end
end
