function OUT = ach_analysis(mat_file, scores_csv, varargin)
% ach_sleep_analysis
% -------------------------------------------------------------
% Analyse ACh (ΔF/F) trace with focus on:
%   1) Sleep transitions (peri-transition ACh dynamics + peak dF/F)
%   2) ACh power spectral density (PSD) during NREM
%
% INPUTS
%   mat_file   : .mat file with ACh trace and its sampling rate.
%                Expected fields (any of these):
%                   ACh trace : 'ach','ACh','ne','dff'
%                   fs (Hz)   : 'ach_frequency','ne_frequency','fs_ach','Fs_ach'
%
%   scores_csv : CSV with 1-Hz sleep scores, MUST have columns:
%                   time_s  (0,1,2,...)
%                   score   (integer code per second)
%
% OPTIONAL NAME/VALUE PAIRS
%   'codes'      : struct with fields .WK .NREM .REM .MA (defaults 1,4,9,15)
%   'epoch_sec'  : epoch length of scoring (default 1 s)
%   'ma_thresh_sec' : Wake runs <= this are re-labeled MA (default 15 s)
%   'reclassify_short_wake_to_MA' : logical (default true)
%
%   'trans_win'  : peri-transition window [pre post] in seconds (default [-100 100])
%   'base_win'   : baseline window (relative to transition) for offset,
%                  default [-60 -10] s
%   'resp_win'   : response window to search for peak dF/F, default [0 60] s
%
%   'psd_fmax'   : max frequency for ACh PSD (Hz), default 0.15
%   'psd_win_sec': Welch window length (sec) for PSD, default 300 s
%
%   'mouse_id'   : string label
%   'session'    : string label
%   'out_prefix' : base name for export files
%   'out_dir'    : directory to save outputs (default: alongside CSV)
%   'verbose'    : true/false (default true)
%
% OUTPUT (struct OUT)
%   OUT.transitions : struct array, one per transition type
%                     .name
%                     .n_events
%                     .t_rel    (time vector, s)
%                     .traces   (n_time x n_events, baseline-corrected)
%                     .mean     (n_time x 1)
%                     .sem      (n_time x 1)
%                     .peaks    (peak dF/F per event)
%
%   OUT.psd         : struct for NREM ACh PSD
%                     .f       (Hz)
%                     .psd     (power)
%                     .band_power
%                     .peak_freq
%                     .peak_amp
%
%   OUT.files.metrics_csv
%   OUT.files.trans_fig
%   OUT.files.psd_fig
%
% -------------------------------------------------------------

%% -------------------- Parse inputs ---------------------------
p = inputParser;
addRequired(p,'mat_file',@ischar);
addRequired(p,'scores_csv',@ischar);

addParameter(p,'codes',struct('WK',1,'NREM',4,'REM',9,'MA',15),@isstruct);
addParameter(p,'epoch_sec',1,@(x)isscalar(x)&&x>0);
addParameter(p,'ma_thresh_sec',15,@(x)isscalar(x)&&x>0);
addParameter(p,'reclassify_short_wake_to_MA',true,@islogical);

addParameter(p,'trans_win',[-100 100],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'base_win',[-60 -10],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'resp_win',[0 60],@(x)isnumeric(x)&&numel(x)==2);

addParameter(p,'psd_fmin',0.01,@(x)isscalar(x)&&x>=0);
addParameter(p,'psd_fmax',0.15,@(x)isscalar(x)&&x>0);
addParameter(p,'psd_win_sec',300,@(x)isscalar(x)&&x>0);

addParameter(p,'mouse_id','',@ischar);
addParameter(p,'session','',@ischar);
addParameter(p,'out_prefix','',@ischar);
addParameter(p,'out_dir','',@ischar);
addParameter(p,'verbose',true,@islogical);

parse(p, mat_file, scores_csv, varargin{:});
S = p.Results;
C = S.codes;
verbose = S.verbose;

if isempty(S.out_dir)
    S.out_dir = fileparts(scores_csv);
    if isempty(S.out_dir), S.out_dir = pwd; end
end
if ~exist(S.out_dir,'dir'); mkdir(S.out_dir); end

OUT = struct();
OUT.success = false;
OUT.files   = struct();

if verbose
    fprintf('\n=== ACh sleep analysis ===\n');
    fprintf('MAT file    : %s\n', mat_file);
    fprintf('Scores CSV  : %s\n', scores_csv);
    fprintf('Output dir  : %s\n', S.out_dir);
end

%% -------------------- Load scores (1 Hz) ---------------------
T = readtable(scores_csv);
assert(all(ismember({'time_s','score'}, T.Properties.VariableNames)), ...
       'Expected columns time_s and score in %s', scores_csv);

t_sec = double(T.time_s(:));   % 0,1,2,...
code  = double(T.score(:));

n_epoch = numel(code);
dur_scores_s = n_epoch * S.epoch_sec;

if verbose
    fprintf('Scores: %d epochs (%.1f min)\n', n_epoch, dur_scores_s/60);
    fprintf('Unique codes: %s\n', mat2str(unique(code)'));
end

% --- Reclassify short Wake bouts as MA (micro-arousals) ------
if S.reclassify_short_wake_to_MA && isfield(C,'MA')
    if verbose
        fprintf('Reclassifying WK runs <= %d s as MA...\n', S.ma_thresh_sec);
    end
    [st,en] = runs_from_codes(code);
    rc   = code(st);
    rdur = t_sec(en) - t_sec(st) + 1;   % duration in seconds
    is_short_wake = (rc == C.WK) & (rdur <= S.ma_thresh_sec);
    for k = find(is_short_wake).'
        code(st(k):en(k)) = C.MA;
    end
end

%% -------------------- Load ACh trace -------------------------
Smat = load(mat_file);

% Candidate variable names for ACh ΔF/F
cand_ach   = {'ach','ACh','ne','dff','DFF'};
cand_fs_ach= {'ach_frequency','ne_frequency','fs_ach','Fs_ach','imaging_frequency'};

ach = [];
ach_name = '';
for k = 1:numel(cand_ach)
    if isfield(Smat, cand_ach{k})
        ach = Smat.(cand_ach{k});
        ach_name = cand_ach{k};
        break;
    end
end
if isempty(ach)
    error('Could not find ACh variable in %s', mat_file);
end
ach = double(ach(:));

fs_ach = [];
fs_name = '';
for k = 1:numel(cand_fs_ach)
    if isfield(Smat, cand_fs_ach{k})
        fs_ach = Smat.(cand_fs_ach{k});
        fs_name = cand_fs_ach{k};
        break;
    end
end
if isempty(fs_ach)
    error('Could not find ACh sampling rate in %s', mat_file);
end
fs_ach = double(fs_ach);

if verbose
    fprintf('Using ACh variable: %s (n = %d samples)\n', ach_name, numel(ach));
    fprintf('ACh sampling rate : %s = %.3f Hz\n', fs_name, fs_ach);
end

%% -------------------- Align ACh with scores ------------------
% We assume ACh trace and scoring start at the same time.
% Build a continuous time axis for ACh and map each sample
% to the corresponding 1-Hz epoch.

n_ach    = numel(ach);
t_ach    = (0:n_ach-1)' / fs_ach;            % seconds from start
sec_idx  = floor(t_ach / S.epoch_sec) + 1;   % epoch index for each sample
valid    = sec_idx >= 1 & sec_idx <= n_epoch;

ach      = ach(valid);
sec_idx  = sec_idx(valid);
state_vec= code(sec_idx);                    % state per ACh sample

if verbose
    fprintf('ACh duration total    : %.1f min\n', n_ach/fs_ach/60);
    fprintf('Overlap with scores   : %.1f min\n', max(t_ach(valid))/60);
end

%% -------------------- PART 1: Sleep transitions ----------------------
% We define three state ONSETS:
%   1) Wake onset:   any NREM/REM -> WK
%   2) NREM onset:   WK or REM    -> NREM
%   3) REM onset:    NREM         -> REM
%
% For each transition type we:
%   - find all transitions on the 1-Hz score vector
%   - extract peri-transition ACh segments
%   - baseline-correct them
%   - compute peak dF/F in a response window

% Template for each transition result
trans_template = struct( ...
    'name',       '', ...
    'n_events',   0, ...
    't_events',   [], ...
    'idx_epochs', [], ...
    't_rel',      [], ...
    'traces',     [], ...
    'mean',       [], ...
    'sem',        [], ...
    'peaks',      []);

% Initialize empty array with correct fields
OUT.transitions = repmat(trans_template,0,1);

% Specify which transitions we want
trans_specs = struct( ...
    'name',       {'Wake_onset','NREM_onset','REM_onset'}, ...
    'from_codes', { [C.NREM C.REM],        [C.WK C.REM],   C.NREM }, ...
    'to_code',    { C.WK,                  C.NREM,         C.REM });

for it = 1:numel(trans_specs)
    spec = trans_specs(it);

    % Find transitions of this type on the 1-Hz scoring
    [t_events, idx_epochs] = find_transitions(t_sec, code, ...
                                              spec.from_codes, spec.to_code);
    
    % Discard events too close to the edges to fit the window
    pre  = S.trans_win(1);
    post = S.trans_win(2);
    keep = (t_events + pre  >= 0) & ...
           (t_events + post <= max(t_ach));
    t_events   = t_events(keep);
    idx_epochs = idx_epochs(keep);
    
    if isempty(t_events)
        if verbose
            fprintf('%s: no usable transitions.\n', spec.name);
        end
        TOUT          = trans_template;
        TOUT.name     = spec.name;
        TOUT.n_events = 0;
        OUT.transitions(end+1) = TOUT; 
        continue;
    end
    
    % Extract peri-event ACh traces around each transition
    [t_rel, traces] = extract_peri_traces(ach, fs_ach, t_events, S.trans_win);
    
    % Baseline-correct each trace using base_win (e.g. -60 to -10 s)
    traces_bc = baseline_correct(traces, t_rel, S.base_win);
    
    % Compute peak dF/F in the response window (e.g. 0–60 s)
    peaks = compute_peaks(traces_bc, t_rel, S.resp_win);
    
    % Mean and SEM across events
    mean_tr = mean(traces_bc, 2);
    sem_tr  = std(traces_bc, 0, 2) / sqrt(size(traces_bc,2));
    
    % Store in OUT
    TOUT               = trans_template;
    TOUT.name          = spec.name;
    TOUT.n_events      = numel(t_events);
    TOUT.t_events      = t_events;
    TOUT.idx_epochs    = idx_epochs;
    TOUT.t_rel         = t_rel;
    TOUT.traces        = traces_bc;
    TOUT.mean          = mean_tr;
    TOUT.sem           = sem_tr;
    TOUT.peaks         = peaks;
    
    OUT.transitions(end+1) = TOUT; 
    
    if verbose
        fprintf('%s: %d transitions kept.\n', spec.name, numel(t_events));
    end
end

%% -------- Plot transition analysis (similar to Fig 3.10) --------------
% Layout: one row per transition type (Wake_onset, NREM_onset, REM_onset)
% Columns:
%   Left  = mean ± SEM peri-event ACh (baseline normalised)
%   Right = peak dF/F for that transition type (points + mean±SEM)

trans_fig = figure('Name','ACh transitions', ...
                   'Color','w', ...
                   'Visible','on');

nT = numel(OUT.transitions);
if nT == 0
    % No transitions at all
    subplot(1,1,1);
    text(0.5,0.5,'No transitions found','HorizontalAlignment','center');
    axis off;
else
    for it = 1:nT
        TOUT = OUT.transitions(it);
        r    = it;   % row index

        % ---------- LEFT: peri-event mean ± SEM ----------
        subplot(nT,2,2*r-1);

        if TOUT.n_events == 0 || isempty(TOUT.t_rel)
            text(0.5,0.5,[TOUT.name ' (no events)'], ...
                 'HorizontalAlignment','center');
            axis off;
        else
            m  = TOUT.mean;
            se = TOUT.sem;
            t  = TOUT.t_rel(:);

            % shaded SEM band
            xx = [t; flipud(t)];
            yy = [m-se; flipud(m+se)];
            fill(xx, yy, [0.8 0.8 0.8], 'EdgeColor','none'); hold on;

            % mean trace
            plot(t, m, 'k','LineWidth',2);

            yline(0,'k:');
            xline(0,'k--');

            xlabel('Time (s)');
            ylabel('ACh dF/F (baseline norm.)');
            title(strrep(TOUT.name,'_',' '));
            xlim(S.trans_win);
            box off;
        end

        % ---------- RIGHT: peak dF/F ----------
        subplot(nT,2,2*r);

        if TOUT.n_events == 0 || isempty(TOUT.peaks)
            text(0.5,0.5,'no peaks','HorizontalAlignment','center');
            axis off;
        else
            % jittered scatter of individual peaks
            xj = 1 + 0.1*(rand(size(TOUT.peaks))-0.5);
            scatter(xj, TOUT.peaks, 35, 'k','filled'); hold on;

            % mean ± SEM
            mu  = mean(TOUT.peaks);
            sem = std(TOUT.peaks)/sqrt(numel(TOUT.peaks));
            errorbar(1, mu, sem, 'k','LineWidth',2,'CapSize',10);

            xlim([0.5 1.5]);
            set(gca,'XTick',1,'XTickLabel',{'peak'});
            ylabel('Peak ACh dF/F');
            title([strrep(TOUT.name,'_',' ') ' peaks']);
            box off;
        end
    end

    sgtitle(sprintf('%s – %s : ACh transitions', S.mouse_id, S.session));
end

% -------- Save transition figure and register in OUT --------
if isempty(S.out_prefix)
    base_trans = sprintf('AChTransitions_%s_%s', ...
                         safe_str(S.mouse_id), safe_str(S.session));
else
    base_trans = sprintf('%s_AChTransitions', safe_str(S.out_prefix));
end
trans_fig_file = fullfile(S.out_dir, [base_trans '.png']);
saveas(trans_fig, trans_fig_file);
OUT.files.trans_fig = trans_fig_file;

%% -------------------- PART 2: ACh PSD in WK / NREM / REM ------------
% We compute PSD of the slow ACh signal separately for each state:
%   - detrend to remove DC
%   - Welch PSD
%   - keep frequencies between psd_fmin and psd_fmax (e.g. 0.01–0.15 Hz)
%   - summary metrics: total power, peak freq, peak amplitude

state_names  = {'Wake','NREM','REM'};
state_codes  = [C.WK,  C.NREM, C.REM];

OUT.psd = struct();

for i = 1:numel(state_names)
    nm   = state_names{i};
    code_i = state_codes(i);

    mask = (state_vec == code_i);
    sig_state = ach(mask);
    Ls = numel(sig_state);

    if Ls == 0
        if verbose
            fprintf('%s ACh: no samples, skipping PSD.\n', nm);
        end
        PSD = struct('f',[],'psd',[],'band_power',NaN,'peak_freq',NaN,'peak_amp',NaN);
    else
        if verbose
            fprintf('%s ACh: %d samples (%.1f min)\n', nm, Ls, Ls/fs_ach/60);
        end

        [f_psd, psd_raw, band_power, peak_f, peak_amp] = ...
            ach_psd_state(sig_state, fs_ach, S.psd_fmin, S.psd_fmax, ...
                          S.psd_win_sec, verbose, ['ACh PSD ' nm]);

        PSD = struct('f',f_psd,'psd',psd_raw, ...
                     'band_power',band_power, ...
                     'peak_freq',peak_f, ...
                     'peak_amp',peak_amp);
    end

    OUT.psd.(nm) = PSD;
end

% -------- Plot PSDs + separate boxplots for metrics -------------------
psd_fig = figure('Name','ACh PSD per state','Color','w','Visible','on');

% Row 1: PSD curves for each state
for i = 1:numel(state_names)
    nm = state_names{i};
    PSD = OUT.psd.(nm);

    subplot(2,3,i);
    if isempty(PSD.f)
        text(0.5,0.5,['No ' nm ' data'], ...
            'HorizontalAlignment','center');
        axis off;
    else
        % normalize to max for plotting (A.U.) but keep raw in OUT
        psd_norm = PSD.psd / max(PSD.psd);
        plot(PSD.f, psd_norm, 'k','LineWidth',1.5);
        xlabel('Frequency (Hz)');
        ylabel('ACh power (A.U.)');
        xlim([S.psd_fmin S.psd_fmax]);
        ylim([0 1.05]);
        title(['ACh PSD – ' nm]);
        box off;
    end
end

% Row 2: three separate boxplots: power, peak freq, peak amp
powers = [OUT.psd.Wake.band_power, OUT.psd.NREM.band_power, OUT.psd.REM.band_power];
pfreqs = [OUT.psd.Wake.peak_freq,  OUT.psd.NREM.peak_freq,  OUT.psd.REM.peak_freq];
pamps  = [OUT.psd.Wake.peak_amp,   OUT.psd.NREM.peak_amp,   OUT.psd.REM.peak_amp];

cats = categorical({'Wake','NREM','REM'});

% (4) Power
subplot(2,3,4);
boxchart(cats', powers');
ylabel('Power (raw units)');
title('ACh PSD power by state');
box off;

% (5) Peak frequency
subplot(2,3,5);
boxchart(cats', pfreqs');
ylabel('Peak frequency (Hz)');
title('ACh PSD peak freq by state');
box off;

% (6) Peak amplitude
subplot(2,3,6);
boxchart(cats', pamps');
ylabel('Peak amplitude (raw units)');
title('ACh PSD peak amp by state');
box off;

sgtitle(sprintf('%s – %s : ACh PSD (Wake / NREM / REM)', S.mouse_id, S.session));

% Save PSD figure
if isempty(S.out_prefix)
    base_psd = sprintf('AChPSD_%s_%s', safe_str(S.mouse_id), safe_str(S.session));
else
    base_psd = sprintf('%s_AChPSD', safe_str(S.out_prefix));
end
psd_fig_file = fullfile(S.out_dir, [base_psd '.png']);
saveas(psd_fig, psd_fig_file);
OUT.files.psd_fig = psd_fig_file;


%% -------------------- Save metrics to CSV --------------------
if isempty(S.out_prefix)
    base_csv = sprintf('AChMetrics_%s_%s', safe_str(S.mouse_id), safe_str(S.session));
else
    base_csv = sprintf('%s_AChMetrics', safe_str(S.out_prefix));
end
csv_file = fullfile(S.out_dir, [base_csv '.csv']);

% -------- Collect per-transition peak summaries into a table ----------
nT = numel(OUT.transitions);   % number of transition types (Wake_onset, NREM_onset, REM_onset)

if nT == 0
    % No transitions at all → make an empty table
    tbl = table();
else
    % Basic columns for each transition type
    trans_names = {OUT.transitions.name}';
    n_events    = arrayfun(@(x)x.n_events, OUT.transitions)';
    peak_means  = arrayfun(@(x) iff(x.n_events>0 && ~isempty(x.peaks), ...
                                    mean(x.peaks), NaN), OUT.transitions)';
    peak_sems   = arrayfun(@(x) iff(x.n_events>0 && ~isempty(x.peaks), ...
                                    std(x.peaks)/sqrt(x.n_events), NaN), OUT.transitions)';

    % Repeat mouse / session for each row
    mouse_col = repmat({S.mouse_id}, nT, 1);
    sess_col  = repmat({S.session},  nT, 1);

    tbl = table(mouse_col, sess_col, ...
                trans_names, n_events, peak_means, peak_sems, ...
        'VariableNames', {'mouse','session','transition','n_events','peak_mean','peak_sem'});

    % Add NREM PSD metrics as extra columns (store them in the first row)
    tbl.ACh_NREM_PSD_power    = NaN(nT,1);
    tbl.ACh_NREM_PSD_peakfreq = NaN(nT,1);
    tbl.ACh_NREM_PSD_peakamp  = NaN(nT,1);

    if isfield(OUT,'psd') && isfield(OUT.psd,'NREM') && ~isempty(OUT.psd.NREM) ...
            && ~isempty(OUT.psd.NREM.f)
        tbl.ACh_NREM_PSD_power(1)    = OUT.psd.NREM.band_power;
        tbl.ACh_NREM_PSD_peakfreq(1) = OUT.psd.NREM.peak_freq;
        tbl.ACh_NREM_PSD_peakamp(1)  = OUT.psd.NREM.peak_amp;
    end
end



writetable(tbl, csv_file);
OUT.files.metrics_csv = csv_file;

if verbose
    fprintf('\nSaved files:\n');

    if isfield(OUT.files,'trans_fig') && ~isempty(OUT.files.trans_fig)
        fprintf('  Transition figure : %s\n', OUT.files.trans_fig);
    else
        fprintf('  Transition figure : (none)\n');
    end

    if isfield(OUT.files,'psd_fig') && ~isempty(OUT.files.psd_fig)
        fprintf('  PSD figure        : %s\n', OUT.files.psd_fig);
    else
        fprintf('  PSD figure        : (none)\n');
    end

    if isfield(OUT.files,'metrics_csv') && ~isempty(OUT.files.metrics_csv)
        fprintf('  Metrics CSV       : %s\n', OUT.files.metrics_csv);
    else
        fprintf('  Metrics CSV       : (none)\n');
    end

    fprintf('Done.\n\n');
end


OUT.success = true;
end 
