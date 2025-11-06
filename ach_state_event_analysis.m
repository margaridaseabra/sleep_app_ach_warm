function ach_state_event_analysis(mat_file, varargin)
% ach_state_event_analysis(mat_file, 'events_csv','events.csv', 'epoch_sec',1)
% mat_file: scored .mat you already use (contains ACh + sleep_scores)
% Optional:
%   'events_csv' : CSV file with columns: label,time_s  (header required)
%   'epoch_sec'  : epoch length used to score sleep (default 1 s)
%   'trend_span' : LOWESS smoothing span in seconds (default 120 s)
%   'win'        : peri-event window [pre post] in seconds (default [-300 600])
%   'baseline'   : baseline window relative to event in seconds (default [-180 -30])
%
% Figures produced:
%   1) Global ACh with LOWESS trend + colored states + event markers
%   2) State-conditioned ACh distributions and means
%   3) Peri-event averages per label (±SEM + baseline-corrected)
%   4) Peri-event split by state at event time
%
% Output tables (in ./_ach_outputs/<basename>):
%   - state_stats.csv           (mean, median, IQR, n)
%   - event_grand_averages.csv  (per label timecourse)
%   - event_peaks.csv           (per event instance peak/area metrics)
%
% Notes:
%   - If your ACh signal needs demod/405-correction, do that upstream.
%   - Sleep state coding expected as integers: WAKE, NREM, REM (map below).
%   - You can hardcode events in the block "INLINE EVENTS" if no CSV.

%% ----------------------- Parse args -----------------------
p = inputParser;
addRequired(p,'mat_file',@ischar);
addParameter(p,'events_csv','',@ischar);
addParameter(p,'epoch_sec',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'trend_span',120,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'win',[-300 600],@(x)isnumeric(x)&&numel(x)==2&&x(1)<x(2));
addParameter(p,'baseline',[-180 -30],@(x)isnumeric(x)&&numel(x)==2&&x(1)<x(2));
parse(p,mat_file,varargin{:});
S = p.Results;

[~,base,~] = fileparts(S.mat_file);
out_dir = fullfile(pwd,'_ach_outputs',base);
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% ----------------------- Load data ------------------------
% Expect variables (rename here if needed):
% ACh vector:       ne OR ach (dFF)             (use one)
% ACh sampling:     ne_frequency OR ach_frequency (Hz)
% Sleep scores vec: sleep_scores (int per epoch)
% EEG sampling:     eeg_frequency (optional)
D = load(S.mat_file);

% --- Resolve ACh channel ---
if isfield(D,'ne'), ach = D.ne(:); 
elseif isfield(D,'ach'), ach = D.ach(:);
else, error('No ACh trace found (ne/ach).');
end

% --- Resolve ACh fs ---
if isfield(D,'ne_frequency'), fs_ach = D.ne_frequency;
elseif isfield(D,'ach_frequency'), fs_ach = D.ach_frequency;
else, error('No ACh sampling frequency (ne_frequency/ach_frequency).');
end
t_ach = (0:numel(ach)-1)'/fs_ach;  % seconds

% --- Sleep labels ---
if ~isfield(D,'sleep_scores'), error('sleep_scores not found.'); end
sleep_scores = D.sleep_scores(:);  % integer codes per epoch
EPOCH_SEC = S.epoch_sec;          % 1s epochs recommended for FP alignment

% Map your codes here if needed:
% Example: 0=Wake, 1=NREM, 2=REM (edit to match your scoring)
WAKE_CODE = 0; NREM_CODE = 1; REM_CODE = 2;

% Build per-sample state masks matched to ACh samples
nS = numel(ach);
nEpoch = numel(sleep_scores);
T_total = nEpoch * EPOCH_SEC;
if abs(T_total - t_ach(end)) > max(1, 0.01*T_total)
    warning('ACh and scoring duration differ (%.1fs vs %.1fs). Truncating to overlap.', t_ach(end), T_total);
end
T = min(t_ach(end), T_total);
idx_keep = t_ach <= T;
ach = ach(idx_keep); t_ach = t_ach(idx_keep); nS = numel(ach);
% epoch index for each ACh sample (1-based)
epoch_idx = min(floor(t_ach/EPOCH_SEC)+1, nEpoch);
state_vec = sleep_scores(epoch_idx);
isWAKE = state_vec==WAKE_CODE;
isNREM = state_vec==NREM_CODE;
isREM  = state_vec==REM_CODE;


%% ----------------- Load/define events ---------------------
% Option A: CSV with header: label,time_s
events = struct('label',{},'time_s',{});
if ~isempty(S.events_csv) && exist(S.events_csv,'file')
    Tcsv = readtable(S.events_csv,'TextType','string');
    assert(all(ismember(["label","time_s"], string(Tcsv.Properties.VariableNames))), ...
        'CSV must have columns: label,time_s');
    for k=1:height(Tcsv)
        if Tcsv.time_s(k) >= 0 && Tcsv.time_s(k) <= T
            events(end+1).label = string(Tcsv.label(k));
            events(end).time_s = double(Tcsv.time_s(k));
        end
    end
else
    % ---------- INLINE EVENTS (edit this block) ----------
    % Example events (seconds since recording start):
    % events = [ struct('label',"ambient_warming_start",'time_s', 3600), ...
    %            struct('label',"saline_injection",'time_s',        5400), ...
    %            struct('label',"drug_injection",'time_s',          7200) ];
    events = events; % leave empty if none
end

% normalize labels to categorical
for k=1:numel(events)
    events(k).label = string(events(k).label);
end
labels = unique(string({events.label}))';

%% ----------------- General trend (LOWESS) -----------------
% Robust LOWESS over time; span in seconds
span_pts = max(5, round(S.trend_span * fs_ach));
trend = smooth(ach, span_pts, 'rlowess');
resid = ach - trend;

%% ----------------- Bout extraction + per-bout stats -----------------
% Epoch-level run-length encode (RLE) of sleep_scores
s = sleep_scores(:);
edges = [true; diff(s)~=0; true];
start_ep = find(edges(1:end-1));
end_ep   = find(edges(2:end)-1);
bout_code = s(start_ep);
nB = numel(start_ep);

% Helper for IQR that ignores NaN
nan_iqr = @(x) diff(quantile(x(~isnan(x)), [0.25 0.75]));

% Prealloc
bout_rows = cell(nB,1);

for k = 1:nB
    code = bout_code(k);
    t0 = (start_ep(k)-1)*EPOCH_SEC;  % seconds (inclusive)
    t1 =  end_ep(k)*EPOCH_SEC;       % seconds (exclusive of next epoch)
    % Map to ACh sample indices
    s0 = max(1, floor(t0*fs_ach)+1);
    s1 = min(nS, max(s0, floor(t1*fs_ach)));

    xi = ach(s0:s1);
    ti = trend(s0:s1);
    tt = t_ach(s0:s1);

    % Robust slope of trend (dF/F per minute)
    if numel(tt) >= 3
        p  = polyfit(tt, ti, 1);       % slope in dF/F per sec
        slope_trend_per_min = p(1)*60; % convert to per-minute
    else
        slope_trend_per_min = NaN;
    end

    % Events inside this bout
    in_labels = strings(0,1);
    for e = 1:numel(events)
        if events(e).time_s >= t0 && events(e).time_s < t1
            in_labels(end+1,1) = events(e).label; %#ok<AGROW>
        end
    end
    if isempty(in_labels), event_labels_in = ""; else, event_labels_in = strjoin(in_labels,';'); end
    n_events_in = numel(in_labels);

    % Build row
    bout_rows{k} = table( ...
        k, code, start_ep(k), end_ep(k), t0, t1, (t1 - t0), ...
        mean(xi,'omitnan'), median(xi,'omitnan'), nan_iqr(xi), std(xi,0,'omitnan'), numel(xi), ...
        mean(ti,'omitnan'), median(ti,'omitnan'), nan_iqr(ti), std(ti,0,'omitnan'), numel(ti), ...
        slope_trend_per_min, string(event_labels_in), n_events_in, ...
        'VariableNames', {'bout_id','state_code','start_epoch','end_epoch','start_time_s','end_time_s','dur_s', ...
                          'ach_mean','ach_median','ach_iqr','ach_std','ach_n', ...
                          'trend_mean','trend_median','trend_iqr','trend_std','trend_n', ...
                          'trend_slope_per_min','event_labels_in','n_events_in'});
end

bout_stats = vertcat(bout_rows{:});

% Optional: human-readable state name
state_name = strings(height(bout_stats),1);
state_name(bout_stats.state_code==WAKE_CODE) = "WAKE";
state_name(bout_stats.state_code==NREM_CODE) = "NREM";
state_name(bout_stats.state_code==REM_CODE ) = "REM";
bout_stats.state = state_name;

% Reorder columns
bout_stats = movevars(bout_stats, 'state', 'After', 'state_code');

% Save
writetable(bout_stats, fullfile(out_dir,'bout_stats.csv'));

%% ----------------- State-conditioned stats ----------------
% Safer IQR for vectors with NaNs
nan_iqr = @(x) diff(quantile(x(~isnan(x)), [0.25 0.75])) ;

mk = {isWAKE, isNREM, isREM};
nm = ["WAKE","NREM","REM"];
rows = cell(3,1);

for i = 1:3
    xi = ach(mk{i});
    ti = trend(mk{i});

    rows{i} = table( ...
        nm(i), ...
        mean(xi,'omitnan'), median(xi,'omitnan'), nan_iqr(xi), numel(xi), ...
        mean(ti,'omitnan'), median(ti,'omitnan'), nan_iqr(ti), numel(ti), ...
        'VariableNames', {'state', ...
                          'ach_mean','ach_median','ach_iqr','ach_n', ...
                          'trend_mean','trend_median','trend_iqr','trend_n'});
end

tab = vertcat(rows{:});
writetable(tab, fullfile(out_dir,'state_stats.csv'));


%% ----------------- Peri-event analysis --------------------
% Window and baseline
tvec = (S.win(1):1:S.win(2))';
nT = numel(tvec);
event_tbl_rows = {};
grand = []; % timecourse per label
peaks_rows = {};

for L = 1:numel(labels)
    lab = labels(L);
    ev_idx = find(strcmp(string({events.label}), lab));
    X = nan(nT, numel(ev_idx));
    X_state = strings(numel(ev_idx),1);
    t_evt = zeros(numel(ev_idx),1);

    for j = 1:numel(ev_idx)
        t0 = events(ev_idx(j)).time_s;
        t_evt(j) = t0;
        idx = round((t0 + tvec)*fs_ach);
        ok = idx >= 1 & idx <= nS;
        x = nan(nT,1); x(ok) = ach(idx(ok));
        % baseline correction
        bmask = tvec>=S.baseline(1) & tvec<=S.baseline(2);
        b = mean(x(bmask),'omitnan');
        X(:,j) = x - b;

        % state at event time (epoch-aligned)
        ep = min(max(1, floor(t0/EPOCH_SEC)+1), nEpoch);
        sc = sleep_scores(ep);
        if sc==WAKE_CODE, X_state(j)="WAKE";
        elseif sc==NREM_CODE, X_state(j)="NREM";
        elseif sc==REM_CODE,  X_state(j)="REM";
        else, X_state(j)="UNK";
        end

        % simple per-event metrics
        post_mask = tvec>=0 & tvec<=300; % first 5 min after event
        [peak_val, peak_idx] = max(x(post_mask),[],'omitnan');
        t_peak = tvec(find(post_mask,1,'first') + peak_idx - 1);
        auc = trapz(tvec(post_mask), x(post_mask));
        peaks_rows(end+1,:) = {char(lab), t0, char(X_state(j)), peak_val, t_peak, auc};
    end

    m = mean(X,2,'omitnan');
    s = std(X,0,2,'omitnan');
    n = sum(~isnan(X),2);
    se = s ./ max(1,sqrt(n));

    % store grand average for export
    grand = [grand; table(repmat(lab,nT,1), tvec, m, se, n, 'VariableNames', ...
        {'label','time_rel_s','mean_dFF_bc','se','n'})];

    % keep all individual per-event timecourses (optional)
    for j=1:numel(ev_idx)
        event_tbl_rows{end+1,1} = table(repmat(lab,nT,1), repmat(t_evt(j),nT,1), tvec, X(:,j), ...
            'VariableNames', {'label','event_time_s','time_rel_s','dFF_bc'});
    end

    % --------- PLOT: peri-event per label ---------
    figure('Color','w','Name',sprintf('Peri-Event ACh — %s', lab));
    plot(tvec, m, 'LineWidth',1.8); hold on;
    fill_between(tvec, m-se, m+se, 0.2);
    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from event (s)'); ylabel('\DeltaF/F (baseline-corr.)');
    title(sprintf('ACh peri-event average — %s', lab));
    grid on; box on;
end

% Export tables
if ~isempty(event_tbl_rows)
    all_ev_tc = vertcat(event_tbl_rows{:});
    writetable(all_ev_tc, fullfile(out_dir,'event_all_timecourses.csv'));
end
if ~isempty(grand)
    writetable(grand, fullfile(out_dir,'event_grand_averages.csv'));
end
if ~isempty(peaks_rows)
    peaks = cell2table(peaks_rows, 'VariableNames', ...
        {'label','event_time_s','state_at_event','peak_val','t_peak_s','auc_0to300s'});
    writetable(peaks, fullfile(out_dir,'event_peaks.csv'));
end

%% ------------ FIG 1 (enhanced): ACh + trend + events + hypnogram -------
figure('Color','w','Name','Whole recording — ACh, states, events');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% --- Top: ACh with LOWESS trend, state bands, event labels ---
ax1 = nexttile;
plot(t_ach, ach, 'LineWidth', 0.6); hold on;
plot(t_ach, trend, 'LineWidth', 1.8);
ylabel('\DeltaF/F'); title('ACh (raw) + LOWESS trend');
grid on; box on;

% Transparent state patches
yl = ylim(ax1);
add_state_patch(t_ach, isWAKE, [0.85 0.92 1.00], yl, 'WAKE');   % pale blue
add_state_patch(t_ach, isNREM, [1.00 0.88 0.88], yl, 'NREM');   % pale red
add_state_patch(t_ach, isREM,  [0.88 1.00 0.88], yl, 'REM');    % pale green

% Event markers + text labels (on the top)
for k=1:numel(events)
    xline(ax1, events(k).time_s, 'k-', 'LineWidth', 1.0);
    text(events(k).time_s, yl(2), sprintf('  %s', events(k).label), ...
        'Parent', ax1,'VerticalAlignment','top','Rotation',90,'FontSize',8);
end
legend(ax1, {'ACh','LOWESS trend'}, 'Location','best');

% --- Bottom: Hypnogram (epoch-level) ---
ax2 = nexttile; hold(ax2,'on');
% Build a step-wise hypnogram aligned to epoch edges
te = (0:nEpoch)*EPOCH_SEC;
ss = [sleep_scores(:); sleep_scores(end)];  % repeat last for stairs
stairs(ax2, te, double(ss), 'LineWidth', 1.2);
ylim(ax2, [min(ss)-0.5, max(ss)+0.5]);
yticks(ax2, unique([WAKE_CODE NREM_CODE REM_CODE]));
yticklabels(ax2, {'WAKE','NREM','REM'});    % order will match your codes
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'State');
title(ax2, 'Hypnogram (epoch scoring)');
grid(ax2,'on'); box(ax2,'on');

linkaxes([ax1 ax2],'x');

%% ------------ FIG 2: State-conditioned distributions --------------------
figure('Color','w','Name','State-conditioned ACh');
subplot(1,2,1);
boxchart(categorical([repmat("WAKE",sum(isWAKE),1); repmat("NREM",sum(isNREM),1); repmat("REM",sum(isREM),1)]), ...
         [ach(isWAKE); ach(isNREM); ach(isREM)]);
ylabel('ACh (dF/F)'); title('Raw ACh by state'); grid on;

subplot(1,2,2);
boxchart(categorical([repmat("WAKE",sum(isWAKE),1); repmat("NREM",sum(isNREM),1); repmat("REM",sum(isREM),1)]), ...
         [trend(isWAKE); trend(isNREM); trend(isREM)]);
ylabel('LOWESS trend (dF/F)'); title('Trend by state'); grid on;

%% ------------ FIG 3: Peri-event split by state at event ------------------
if ~isempty(events)
    figure('Color','w','Name','Peri-event by state-at-event');
    states_at_evt = ["WAKE","NREM","REM"];
    for s=1:3
        st = states_at_evt(s);
        hold on;
        rows = [];
        for L=1:numel(labels)
            lab = labels(L);
            sub = grand(strcmp(grand.label, lab),:); %#ok<NODEF>
            % For legend clarity, overlay each label as separate line
            plot(sub.time_rel_s, sub.mean_dFF_bc, 'LineWidth', 1.2, 'DisplayName', sprintf('%s (all)', lab));
        end
    end
    xline(0,'k-'); yline(0,'k:'); grid on; box on;
    xlabel('Time from event (s)'); ylabel('\DeltaF/F (baseline-corr.)');
    title('Peri-event (all labels combined): mean \pm SE');
    legend('Location','bestoutside');
end

%% ---------------------- Console summary --------------------
fprintf('\nSaved outputs in: %s\n', out_dir);
disp(tab);

end % main fn

% ===================== helpers ======================
function fill_between(x, y1, y2, alpha)
    if any(isnan(y1)) || any(isnan(y2)), hold on; return; end
    xv = [x; flipud(x)];
    yv = [y1; flipud(y2)];
    p = patch(xv, yv, [0 0 0], 'EdgeColor','none');
    if verLessThan('matlab','9.5'), set(p,'FaceAlpha',alpha); else, p.FaceAlpha = alpha; end
end

function add_state_patch(t, mask, color, yl, ~)
    dd = diff([false; mask; false]);
    on  = find(dd==1);
    off = find(dd==-1)-1;
    for i=1:numel(on)
        x0 = t(on(i)); x1 = t(off(i));
        patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], color, ...
              'EdgeColor','none', 'FaceAlpha', 0.12);
    end
end
