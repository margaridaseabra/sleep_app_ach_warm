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
%   'slope_win'  : window for transition slope fit (default [0 60] s)
%
%   'psd_fmin'   : min frequency for ACh PSD (Hz), default 0.01
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
%                     .t_events
%                     .idx_epochs
%                     .t_rel    (time vector, s)
%                     .traces   (n_time x n_events, baseline-corrected)
%                     .mean     (n_time x 1)
%                     .sem      (n_time x 1)
%                     .peaks    (peak dF/F per event)
%                     .slopes   (slope per event in slope_win)
%
%   OUT.psd         : struct with fields .Wake .NREM .REM
%                     each: .f, .psd, .band_power, .peak_freq, .peak_amp
%
%   OUT.state_slopes: struct with fields .Wake .NREM .REM
%                     each: .mean, .sem, .n
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
addParameter(p,'slope_win',[0 60],@(x)isnumeric(x)&&numel(x)==2);

addParameter(p,'psd_fmin',0.01,@(x)isscalar(x)&&x>=0);
addParameter(p,'psd_fmax',0.15,@(x)isscalar(x)&&x>0);
addParameter(p,'psd_win_sec',300,@(x)isscalar(x)&&x>0);

addParameter(p,'mouse_id','',@ischar);
addParameter(p,'session','',@ischar);
addParameter(p,'out_prefix','',@ischar);
addParameter(p,'out_dir','',@ischar);
addParameter(p,'verbose',true,@islogical);

parse(p, mat_file, scores_csv, varargin{:});
S        = p.Results;
C        = S.codes;
verbose  = S.verbose;
slope_win = S.slope_win;

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

n_epoch      = numel(code);
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
cand_ach    = {'ach','ACh','ne','dff','DFF'};
cand_fs_ach = {'ach_frequency','ne_frequency','fs_ach','Fs_ach','imaging_frequency'};

ach      = [];
ach_name = '';
for k = 1:numel(cand_ach)
    if isfield(Smat, cand_ach{k})
        ach      = Smat.(cand_ach{k});
        ach_name = cand_ach{k};
        break;
    end
end
if isempty(ach)
    error('Could not find ACh variable in %s', mat_file);
end
ach = double(ach(:));

fs_ach   = [];
fs_name  = '';
for k = 1:numel(cand_fs_ach)
    if isfield(Smat, cand_fs_ach{k})
        fs_ach  = Smat.(cand_fs_ach{k});
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
n_ach    = numel(ach);
t_ach    = (0:n_ach-1)' / fs_ach;            % seconds from start
sec_idx  = floor(t_ach / S.epoch_sec) + 1;   % epoch index for each sample
valid    = sec_idx >= 1 & sec_idx <= n_epoch;

ach       = ach(valid);
sec_idx   = sec_idx(valid);
state_vec = code(sec_idx);                   % state per ACh sample

if verbose
    fprintf('ACh duration total    : %.1f min\n', n_ach/fs_ach/60);
    fprintf('Overlap with scores   : %.1f min\n', max(t_ach(valid))/60);
end

%% -------------------- PART 1: Sleep transitions ----------------------
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
    'peaks',      [], ...
    'slopes',     [] );

% Initialize empty array with correct fields
OUT.transitions = repmat(trans_template,0,1);

% Specify which transitions we want
trans_specs = struct( ...
    'name',       {'Wake_onset','NREM_onset','REM_onset'}, ...
    'from_codes', { [C.NREM C.REM],         [C.WK C.REM],  C.NREM }, ...
    'to_code',    { C.WK,                   C.NREM,        C.REM });

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

    % ----------- NON-EMPTY CASE: compute everything ----------------
    % Extract peri-event ACh traces around each transition
    [t_rel, traces] = extract_peri_traces(ach, fs_ach, t_events, S.trans_win);
    
    % Baseline-correct each trace using base_win (e.g. -60 to -10 s)
    traces_bc = baseline_correct(traces, t_rel, S.base_win);
    
    % Compute peak dF/F in the response window (e.g. 0–60 s)
    peaks = compute_peaks(traces_bc, t_rel, S.resp_win);

    % Compute slope per event in slope_win
    slopes = compute_slopes(traces_bc, t_rel, slope_win);
    
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
    TOUT.slopes        = slopes;
    
    OUT.transitions(end+1) = TOUT; 
    
    if verbose
        fprintf('%s: %d transitions kept.\n', spec.name, numel(t_events));
    end
end

%% -------- Plot transition analysis (similar to Fig 3.10) --------------
trans_fig = figure('Name','ACh transitions', ...
                   'Color','w', ...
                   'Visible','on');

nT = numel(OUT.transitions);
if nT == 0
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

            xx = [t; flipud(t)];
            yy = [m-se; flipud(m+se)];
            fill(xx, yy, [0.8 0.8 0.8], 'EdgeColor','none'); hold on;

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
            xj = 1 + 0.1*(rand(size(TOUT.peaks))-0.5);
            scatter(xj, TOUT.peaks, 35, 'k','filled'); hold on;

            mu  = mean(TOUT.peaks);
            se  = std(TOUT.peaks)/sqrt(numel(TOUT.peaks));
            errorbar(1, mu, se, 'k','LineWidth',2,'CapSize',10);

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
state_names = {'Wake','NREM','REM'};
state_codes = [C.WK,  C.NREM, C.REM];

OUT.psd = struct();

for i = 1:numel(state_names)
    nm     = state_names{i};
    code_i = state_codes(i);

    mask      = (state_vec == code_i);
    sig_state = ach(mask);
    Ls        = numel(sig_state);

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

%% -------------------- State-wise ACh slope -------------------
OUT.state_slopes = struct();

dt   = 1/fs_ach;
dach = gradient(ach) / dt;     % numerical derivative (dF/F per s)

for i = 1:numel(state_names)
    nm     = state_names{i};
    code_i = state_codes(i);

    mask      = (state_vec == code_i);
    slopes_i  = dach(mask);
    nMask     = nnz(mask);

    if nMask < 5
        OUT.state_slopes.(nm).mean = NaN;
        OUT.state_slopes.(nm).sem  = NaN;
        OUT.state_slopes.(nm).n    = nMask;
    else
        OUT.state_slopes.(nm).mean = mean(slopes_i,'omitnan');
        OUT.state_slopes.(nm).sem  = std(slopes_i,'omitnan') / sqrt(nMask);
        OUT.state_slopes.(nm).n    = nMask;
    end
end

%% -------- Plot PSDs + simple summary metrics ----------------
psd_fig = figure('Name','ACh PSD per state','Color','w','Visible','on');

% Row 1: PSD curves for each state
for i = 1:numel(state_names)
    nm  = state_names{i};
    PSD = OUT.psd.(nm);

    subplot(2,3,i);
    if isempty(PSD.f)
        text(0.5,0.5,['No ' nm ' data'], ...
            'HorizontalAlignment','center');
        axis off;
    else
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
%% -------------------- Slope plots (transitions + states) -------------
slope_fig = figure('Name','ACh slopes','Color','w','Visible','on');

nT = numel(OUT.transitions);

% ---------- ROW 1: per-transition slopes (scatter + mean ± SEM) ------
for it = 1:nT
    TOUT = OUT.transitions(it);
    subplot(2, max(nT,3), it);  % keep at least 3 columns so layout isn't crazy

    if TOUT.n_events == 0 || isempty(TOUT.slopes)
        text(0.5,0.5,[TOUT.name ' (no events)'], ...
            'HorizontalAlignment','center');
        axis off;
    else
        xj = 1 + 0.1*(rand(size(TOUT.slopes))-0.5);
        scatter(xj, TOUT.slopes, 30, 'k','filled'); hold on;

        mu  = mean(TOUT.slopes, 'omitnan');
        se  = std(TOUT.slopes, 'omitnan') / sqrt(sum(~isnan(TOUT.slopes)));
        errorbar(1, mu, se, 'k','LineWidth',1.5,'CapSize',8);

        xlim([0.5 1.5]);
        set(gca,'XTick',1,'XTickLabel',{'slope'});
        ylabel('Slope (dF/F per s)');
        title(strrep(TOUT.name,'_',' '));
        box off;
    end
end

% ---------- ROW 2: state-wise mean slopes (Wake / NREM / REM) --------
subplot(2, max(nT,3), max(nT,3) + 1 : 2*max(nT,3));  % span entire second row
hold on;

state_names = {'Wake','NREM','REM'};
mu_s = nan(1,numel(state_names));
se_s = nan(1,numel(state_names));

for i = 1:numel(state_names)
    nm = state_names{i};
    if isfield(OUT.state_slopes, nm)
        mu_s(i) = OUT.state_slopes.(nm).mean;
        se_s(i) = OUT.state_slopes.(nm).sem;
    end
end

x = 1:numel(state_names);
bar(x, mu_s, 'FaceColor',[0.6 0.6 0.6]); 
errorbar(x, mu_s, se_s, 'k','LineStyle','none','CapSize',8);

set(gca,'XTick',x,'XTickLabel',state_names);
ylabel('Slope (dF/F per s)');
title('Mean ACh slope by state');
box off;

sgtitle(sprintf('%s – %s : ACh slopes', S.mouse_id, S.session));

% Save slope figure
if isempty(S.out_prefix)
    base_slope = sprintf('AChSlopes_%s_%s', safe_str(S.mouse_id), safe_str(S.session));
else
    base_slope = sprintf('%s_AChSlopes', safe_str(S.out_prefix));
end
slope_fig_file = fullfile(S.out_dir, [base_slope '.png']);
saveas(slope_fig, slope_fig_file);
OUT.files.slope_fig = slope_fig_file;

%% -------------------- Save metrics to CSV --------------------
if isempty(S.out_prefix)
    base_csv = sprintf('AChMetrics_%s_%s', safe_str(S.mouse_id), safe_str(S.session));
else
    base_csv = sprintf('%s_AChMetrics', safe_str(S.out_prefix));
end
csv_file = fullfile(S.out_dir, [base_csv '.csv']);

nT = numel(OUT.transitions);

if nT == 0
    tbl = table();
else
    trans_names = {OUT.transitions.name}';
    n_events    = arrayfun(@(x)x.n_events, OUT.transitions)';
    peak_means  = arrayfun(@(x) iff(x.n_events>0 && ~isempty(x.peaks), ...
                                    mean(x.peaks), NaN), OUT.transitions)';
    peak_sems   = arrayfun(@(x) iff(x.n_events>0 && ~isempty(x.peaks), ...
                                    std(x.peaks)/sqrt(x.n_events), NaN), OUT.transitions)';

    slope_means = arrayfun(@(x) iff(x.n_events>0 && isfield(x,'slopes') ...
                                    && ~isempty(x.slopes), ...
                                    mean(x.slopes,'omitnan'), NaN), ...
                           OUT.transitions)';
    slope_sems  = arrayfun(@(x) iff(x.n_events>0 && isfield(x,'slopes') ...
                                    && ~isempty(x.slopes), ...
                                    std(x.slopes,'omitnan')/sqrt(x.n_events), NaN), ...
                           OUT.transitions)';

    mouse_col = repmat({S.mouse_id}, nT, 1);
    sess_col  = repmat({S.session},  nT, 1);

    tbl = table(mouse_col, sess_col, ...
                trans_names, n_events, ...
                peak_means, peak_sems, ...
                slope_means, slope_sems, ...
       'VariableNames', {'mouse','session','transition','n_events', ...
                         'peak_mean','peak_sem', ...
                         'slope_mean','slope_sem'});

    % NREM PSD metrics in first row
    tbl.ACh_NREM_PSD_power    = NaN(nT,1);
    tbl.ACh_NREM_PSD_peakfreq = NaN(nT,1);
    tbl.ACh_NREM_PSD_peakamp  = NaN(nT,1);

    if isfield(OUT,'psd') && isfield(OUT.psd,'NREM') && ~isempty(OUT.psd.NREM) ...
            && ~isempty(OUT.psd.NREM.f)
        tbl.ACh_NREM_PSD_power(1)    = OUT.psd.NREM.band_power;
        tbl.ACh_NREM_PSD_peakfreq(1) = OUT.psd.NREM.peak_freq;
        tbl.ACh_NREM_PSD_peakamp(1)  = OUT.psd.NREM.peak_amp;
    end

    % state-wise mean ACh slope (Wake, NREM, REM) in first row
    tbl.ACh_slope_Wake = NaN(nT,1);
    tbl.ACh_slope_NREM = NaN(nT,1);
    tbl.ACh_slope_REM  = NaN(nT,1);

    if isfield(OUT,'state_slopes')
        if isfield(OUT.state_slopes,'Wake')
            tbl.ACh_slope_Wake(1) = OUT.state_slopes.Wake.mean;
        end
        if isfield(OUT.state_slopes,'NREM')
            tbl.ACh_slope_NREM(1) = OUT.state_slopes.NREM.mean;
        end
        if isfield(OUT.state_slopes,'REM')
            tbl.ACh_slope_REM(1)  = OUT.state_slopes.REM.mean;
        end
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

%% ======================= Helper: compute_slopes =======================
function slopes = compute_slopes(traces, t_rel, win)
% Compute per-event linear slope of baseline-corrected traces
% within a specified time window WIN = [t_start t_end] (s)
%
% traces: [nTime x nEvents]
% t_rel : [nTime x 1] time vector (s)

if isempty(traces)
    slopes = [];
    return;
end

t    = t_rel(:);
mask = (t >= win(1)) & (t <= win(2));
t_win = t(mask);

nEv    = size(traces,2);
slopes = nan(nEv,1);

for j = 1:nEv
    y = traces(mask,j);
    if all(isnan(y)) || numel(y) < 3
        slopes(j) = NaN;
        continue;
    end
    p = polyfit(t_win, y, 1);     % y = p(1)*t + p(2)
    slopes(j) = p(1);             % slope (dF/F per second)
end
end
