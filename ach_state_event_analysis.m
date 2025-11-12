function ach_state_event_analysis(mat_file, varargin)
% ach_state_event_analysis(mat_file, 'events_csv','events.csv', 'epoch_sec',1)
% mat_file: scored .mat you already use (contains ACh + sleep_scores)
% Optional:
%   'events_csv' : CSV file with columns: label,time_s  (header required)
%   'epoch_sec'  : epoch length used to score sleep (default 1 s)
%   'trend_span' : LOWESS smoothing span in seconds (default 120 s)
%   'win'        : peri-event window [pre post] in seconds (default [-300 600])
%   'baseline'   : baseline window relative to event in seconds (default [-180 -30])

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
fig_dir = fullfile(out_dir,'figs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% ----------------------- Load data ------------------------
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
sleep_scores = D.sleep_scores(:);
EPOCH_SEC = S.epoch_sec;

% Map your codes here:
WAKE_CODE = 0; NREM_CODE = 1; REM_CODE = 2;

% Align scoring with ACh samples
nS = numel(ach);
nEpoch = numel(sleep_scores);
T_total = nEpoch * EPOCH_SEC;
if abs(T_total - t_ach(end)) > max(1, 0.01*T_total)
    warning('ACh and scoring duration differ (%.1fs vs %.1fs). Truncating to overlap.', t_ach(end), T_total);
end
T = min(t_ach(end), T_total);
idx_keep = t_ach <= T;
ach = ach(idx_keep); t_ach = t_ach(idx_keep); nS = numel(ach);
epoch_idx = min(floor(t_ach/EPOCH_SEC)+1, nEpoch);
state_vec = sleep_scores(epoch_idx);
isWAKE = state_vec==WAKE_CODE;
isNREM = state_vec==NREM_CODE;
isREM  = state_vec==REM_CODE;

%% ----------------- Load/define events ---------------------
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
end
for k=1:numel(events), events(k).label = string(events(k).label); end
labels = unique(string({events.label}))';

%% ----------------- NaN-safe helpers -----------------
nan_mean = @(x) mean(x(~isnan(x)));
nan_med  = @(x) median(x(~isnan(x)));
nan_std  = @(x) std(x(~isnan(x)));
nan_iqr  = @(x) diff(quantile(x(~isnan(x)), [0.25 0.75]));
nan_rms  = @(x) sqrt(mean((x(~isnan(x))).^2));
nan_mad  = @(x) nan_med(abs(x - nan_med(x)));
zscore_n = @(x) (x - nan_mean(x)) ./ nan_std(x);

nan_skew = @(x) (numel(x(~isnan(x)))>=3) .* skewness(x(~isnan(x))) + ...
                (numel(x(~isnan(x)))<3)  .* NaN;
nan_kurt = @(x) (numel(x(~isnan(x)))>=4) .* kurtosis(x(~isnan(x))) + ...
                (numel(x(~isnan(x)))<4)  .* NaN;

%% ----------------- General trend (LOWESS) -----------------
span_pts = max(5, round(S.trend_span * fs_ach));
ach_trend = smooth(ach, span_pts, 'rlowess');
resid = ach - ach_trend;

%% ----------------- Bout extraction + per-bout stats -----------------
s = sleep_scores(:);
edges = [true; diff(s)~=0; true];
start_ep = find(edges(1:end-1));
end_ep   = find(edges(2:end)-1);
bout_code = s(start_ep);
nB = numel(start_ep);

bout_rows = cell(nB,1);
for k = 1:nB
    code = bout_code(k);
    t0 = (start_ep(k)-1)*EPOCH_SEC;
    t1 =  end_ep(k)*EPOCH_SEC;
    s0 = max(1, floor(t0*fs_ach)+1);
    s1 = min(nS, max(s0, floor(t1*fs_ach)));

    xi = ach(s0:s1);
    ti = ach_trend(s0:s1);
    tt = t_ach(s0:s1);

    if numel(tt) >= 3
        p  = polyfit(tt, ti, 1);
        slope_ach_trend_per_min = p(1)*60;
    else
        slope_ach_trend_per_min = NaN;
    end

    in_labels = strings(0,1);
    for e = 1:numel(events)
        if events(e).time_s >= t0 && events(e).time_s < t1
            in_labels(end+1,1) = events(e).label; %#ok<AGROW>
        end
    end
    event_labels_in = ""; if ~isempty(in_labels), event_labels_in = strjoin(in_labels,';'); end
    n_events_in = numel(in_labels);

    bout_rows{k} = table( ...
        k, code, start_ep(k), end_ep(k), t0, t1, (t1 - t0), ...
        nan_mean(xi), nan_med(xi), nan_iqr(xi), nan_std(xi), numel(xi), ...
        nan_mean(ti), nan_med(ti), nan_iqr(ti), nan_std(ti), numel(ti), ...
        slope_ach_trend_per_min, string(event_labels_in), n_events_in, ...
        'VariableNames', {'bout_id','state_code','start_epoch','end_epoch','start_time_s','end_time_s','dur_s', ...
                          'ach_mean','ach_median','ach_iqr','ach_std','ach_n', ...
                          'ach_trend_mean','ach_trend_median','ach_trend_iqr','ach_trend_std','ach_trend_n', ...
                          'ach_trend_slope_per_min','event_labels_in','n_events_in'});
end

bout_stats = vertcat(bout_rows{:});
state_name = strings(height(bout_stats),1);
state_name(bout_stats.state_code==WAKE_CODE) = "WAKE";
state_name(bout_stats.state_code==NREM_CODE) = "NREM";
state_name(bout_stats.state_code==REM_CODE ) = "REM";
bout_stats.state = state_name;
bout_stats = movevars(bout_stats, 'state', 'After', 'state_code');
writetable(bout_stats, fullfile(out_dir,'bout_stats.csv'));

%% ----------------- State-conditioned stats ----------------
mk = {isWAKE, isNREM, isREM};
nm = ["WAKE","NREM","REM"];
rows = cell(3,1);
for i = 1:3
    xi = ach(mk{i});
    ti = ach_trend(mk{i});
    rows{i} = table( ...
        nm(i), ...
        nan_mean(xi), nan_med(xi), nan_iqr(xi), numel(xi), ...
        nan_mean(ti), nan_med(ti), nan_iqr(ti), numel(ti), ...
        'VariableNames', {'state', ...
                          'ach_mean','ach_median','ach_iqr','ach_n', ...
                          'ach_trend_mean','ach_trend_median','ach_trend_iqr','ach_trend_n'});
end
tab = vertcat(rows{:});
writetable(tab, fullfile(out_dir,'state_stats.csv'));

%% ----------------- Peri-event analysis --------------------
tvec = (S.win(1):1:S.win(2))';
nT = numel(tvec);
event_tbl_rows = {};
grand = []; peaks_rows = {};

if ~isempty(events)
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
            bmask = tvec>=S.baseline(1) & tvec<=S.baseline(2);
            b = nan_mean(x(bmask));
            X(:,j) = x - b;

            ep = min(max(1, floor(t0/EPOCH_SEC)+1), nEpoch);
            sc = sleep_scores(ep);
            if sc==WAKE_CODE, X_state(j)="WAKE";
            elseif sc==NREM_CODE, X_state(j)="NREM";
            elseif sc==REM_CODE,  X_state(j)="REM";
            else, X_state(j)="UNK";
            end

            post_mask = tvec>=0 & tvec<=300;
            [peak_val, peak_idx] = max(x(post_mask),[],'omitnan');
            t_peak = tvec(find(post_mask,1,'first') + peak_idx - 1);
            auc = trapz(tvec(post_mask), x(post_mask));
            peaks_rows(end+1,:) = {char(lab), t0, char(X_state(j)), peak_val, t_peak, auc};
        end

        m = mean(X,2,'omitnan');
        s = std(X,0,2,'omitnan');
        n = sum(~isnan(X),2);
        se = s ./ max(1,sqrt(n));

        grand = [grand; table(repmat(lab,nT,1), tvec, m, se, n, ...
            'VariableNames', {'label','time_rel_s','mean_dFF_bc','se','n'})];

        for j=1:numel(ev_idx)
            event_tbl_rows{end+1,1} = table(repmat(lab,nT,1), repmat(t_evt(j),nT,1), tvec, X(:,j), ...
                'VariableNames', {'label','event_time_s','time_rel_s','dFF_bc'});
        end

        % Per-label plot
        figure('Color','w','Name',sprintf('Peri-Event ACh — %s', lab));
        plot(tvec, m, 'LineWidth',1.8); hold on;
        fill_between(tvec, m-se, m+se, 0.2);
        xline(0,'k-'); yline(0,'k:');
        xlabel('Time from event (s)'); ylabel('\DeltaF/F (baseline-corr.)');
        title(sprintf('ACh peri-event average — %s', lab));
        grid on; box on;

    end
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

% Combined grand-average plot (only if we have events)
if ~isempty(grand)
    labs = unique(grand.label);
    figure('Color','w','Name','Peri-event: mean \pm SE by label'); hold on; grid on; box on;
    h = gobjects(numel(labs),1);
    for i = 1:numel(labs)
        sub = grand(grand.label==labs(i), :);
        t  = sub.time_rel_s; m = sub.mean_dFF_bc; se = sub.se;
        fill([t; flipud(t)], [m-se; flipud(m+se)], [0 0 0], ...
             'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
        h(i) = plot(t, m, 'LineWidth',1.8, 'DisplayName', char(labs(i)));
    end
    xline(0,'k-','HandleVisibility','off'); yline(0,'k:','HandleVisibility','off');
    xlabel('Time from event (s)'); ylabel('\DeltaF/F (baseline-corr.)');
    title('Peri-event averages (mean \pm SE)'); legend(h,'Location','best');
    % SAVE:
    save_fig(f, fig_dir, "peri_event_all_labels");
   
end

%% ================= ACh trace–centric statistics =================
% Global descriptors
S_global = table( ...
    nan_mean(ach),  nan_med(ach),  nan_iqr(ach),  nan_std(ach),  nan_mad(ach), nan_rms(ach), ...
    nan_skew(ach),  nan_kurt(ach), ...
    nan_mean(ach_trend), nan_med(ach_trend), nan_iqr(ach_trend), nan_std(ach_trend), ...
    nan_mean(resid),     nan_iqr(resid),     nan_std(resid),     nan_mad(resid), ...
    'VariableNames', { ...
      'ach_mean','ach_median','ach_iqr','ach_std','ach_mad','ach_rms','ach_skew','ach_kurt', ...
      'ach_trend_mean','ach_trend_median','ach_trend_iqr','ach_trend_std', ...
      'resid_mean','resid_iqr','resid_std','resid_mad' ...
});


% Autocorr half-life on residual
maxLagSec = 60; L = min(numel(resid), round(maxLagSec*fs_ach));
r = xcorr(resid - nan_mean(resid), L, 'coeff'); r = r(L+1:end);
halflife_s = NaN;
idx = find(r <= 0.5, 1, 'first');
if ~isempty(idx) && idx>=2
    x1 = (idx-2)/fs_ach; x2 = (idx-1)/fs_ach;
    y1 = r(idx-1);       y2 = r(idx);
    halflife_s = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1);
end

% PSD slope (0.005–0.1 Hz)
nfft = 2^nextpow2(round(60*fs_ach));
[pxx, f] = pwelch(ach, hamming(round(60*fs_ach)), [], nfft, fs_ach, 'onesided');
band = f>=0.005 & f<=0.1 & pxx>0;
psd_slope = NaN; psd_intercept = NaN;
if any(band)
    X = [ones(sum(band),1), log10(f(band))]; y = log10(pxx(band)); b = X\y;
    psd_intercept = b(1); psd_slope = b(2);
end
S_times = table(halflife_s, psd_slope, psd_intercept);

% Drift (global slope of trend)
t_s = (0:numel(ach_trend)-1)'/fs_ach;
if numel(t_s) >= 3
    p = polyfit(t_s, ach_trend, 1);
    slope_ach_trend_per_min_global = p(1)*60;
else
    slope_ach_trend_per_min_global = NaN;
end
S_drift = table(slope_ach_trend_per_min_global);

% Transient detection on residual
z = zscore_n(resid);
minZ   = 2.5;
minSep = round(2*fs_ach);
minW   = round(0.3*fs_ach);
[pk, loc, w, prom] = findpeaks(z, 'MinPeakHeight',minZ, 'MinPeakDistance',minSep, ...
                                  'MinPeakProminence',1, 'MinPeakWidth',minW);
peak_time_s = loc / fs_ach;
peak_val_df = resid(loc);
trans_tbl = table(peak_time_s, pk, prom, w/fs_ach, peak_val_df, ...
    'VariableNames', {'t_s','z_peak','prominence_z','width_s','resid_peak_df'});

% Change-points (try/catch if toolbox missing)
cp_tbl = table();
try
    maxN = min(12, max(2, floor(numel(ach_trend)/(fs_ach*60*5))));
    if maxN >= 2
        idx_cp = findchangepts(ach_trend,'Statistic','mean','MaxNumChanges',maxN);
        cp_tbl = [cp_tbl; table((idx_cp(:)-1)/fs_ach, repmat("mean",numel(idx_cp),1), ...
                 'VariableNames', {'t_s','type'})];
    end
    maxN2 = min(12, max(2, floor(numel(resid)/(fs_ach*60*5))));
    if maxN2 >= 2
        idx_cp2 = findchangepts(resid,'Statistic','variance','MaxNumChanges',maxN2);
        cp_tbl = [cp_tbl; table((idx_cp2(:)-1)/fs_ach, repmat("variance",numel(idx_cp2),1), ...
                 'VariableNames', {'t_s','type'})];
    end
catch
end

% Per-state trace-centric stats
rows = cell(3,1);
for i=1:3
    xi  = ach(mk{i});
    tri = ach_trend(mk{i});
    rei = resid(mk{i});
    ti  = t_s(mk{i});
    slope_i = NaN;
    if numel(ti)>=3, p = polyfit(ti, tri, 1); slope_i = p(1)*60; end
    rows{i} = table(nm(i), ...
        nan_mean(xi), nan_med(xi), nan_iqr(xi), nan_std(xi), ...
        nan_mean(tri), nan_med(tri), nan_iqr(tri), nan_std(tri), slope_i, ...
        nan_std(rei), nan_mad(rei), ...
        'VariableNames', {'state', ...
          'ach_mean','ach_median','ach_iqr','ach_std', ...
          'ach_trend_mean','ach_trend_median','ach_trend_iqr','ach_trend_std','ach_trend_slope_per_min', ...
          'resid_std','resid_mad'});
end
S_state = vertcat(rows{:});

% Save trace-centric outputs
writetable([S_global S_times S_drift], fullfile(out_dir,'ach_trace_stats.csv'));
if ~isempty(trans_tbl), writetable(trans_tbl, fullfile(out_dir,'ach_transients.csv')); end
if ~isempty(cp_tbl), writetable(cp_tbl, fullfile(out_dir,'ach_changepoints.csv')); end
writetable(S_state, fullfile(out_dir,'ach_trace_stats_by_state.csv'));
if any(band)
    Tpsd = table(f(band), pxx(band), 'VariableNames', {'f_Hz','Pxx'});
    writetable(Tpsd, fullfile(out_dir,'ach_psd.csv'));
end

%% ------------ FIG 1 (enhanced): ACh + trend + events + hypnogram -------
figure('Color','w','Name','Whole recording — ACh, states, events');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
save_fig(f, fig_dir, "whole_recording");


ax1 = nexttile;
plot(t_ach, ach, 'LineWidth', 0.6); hold on;
plot(t_ach, ach_trend, 'LineWidth', 1.8);
ylabel('\DeltaF/F'); title('ACh (raw) + LOWESS trend'); grid on; box on;
yl = ylim(ax1);
add_state_patch(t_ach, isWAKE, [0.85 0.92 1.00], yl, 'WAKE');
add_state_patch(t_ach, isNREM, [1.00 0.88 0.88], yl, 'NREM');
add_state_patch(t_ach, isREM,  [0.88 1.00 0.88], yl, 'REM');
for k=1:numel(events)
    xline(ax1, events(k).time_s, 'k-', 'LineWidth', 1.0);
    text(events(k).time_s, yl(2), sprintf('  %s', events(k).label), ...
        'Parent', ax1,'VerticalAlignment','top','Rotation',90,'FontSize',8);
end
legend(ax1, {'ACh','LOWESS trend'}, 'Location','best');

ax2 = nexttile; hold(ax2,'on');
te = (0:nEpoch)*EPOCH_SEC;
ss = [sleep_scores(:); sleep_scores(end)];
stairs(ax2, te, double(ss), 'LineWidth', 1.2);
ylim(ax2, [min(ss)-0.5, max(ss)+0.5]);
yticks(ax2, unique([WAKE_CODE NREM_CODE REM_CODE]));
yticklabels(ax2, {'WAKE','NREM','REM'});
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'State'); title(ax2, 'Hypnogram'); grid(ax2,'on'); box(ax2,'on');
linkaxes([ax1 ax2],'x');

%% ------------ FIG 2: State-conditioned distributions --------------------
figure('Color','w','Name','State-conditioned ACh');
subplot(1,2,1);
boxchart(categorical([repmat("WAKE",sum(isWAKE),1); repmat("NREM",sum(isNREM),1); repmat("REM",sum(isREM),1)]), ...
         [ach(isWAKE); ach(isNREM); ach(isREM)]);
ylabel('ACh (dF/F)'); title('Raw ACh by state'); grid on;

subplot(1,2,2);
boxchart(categorical([repmat("WAKE",sum(isWAKE),1); repmat("NREM",sum(isNREM),1); repmat("REM",sum(isREM),1)]), ...
         [ach_trend(isWAKE); ach_trend(isNREM); ach_trend(isREM)]);
ylabel('LOWESS trend (dF/F)'); title('Trend by state'); grid on;
save_fig(f, fig_dir, "state_conditioned");


%% ---------------------- Console summary --------------------
fprintf('\nSaved outputs in: %s\n', out_dir);
disp(tab);

end % main fn

% ===================== helpers ======================
function fill_between(x, y1, y2, alpha)
    % Shade area between two curves
    x = x(:); y1 = y1(:); y2 = y2(:);
    if any(isnan(y1)) || any(isnan(y2)), return; end
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
function save_fig(fig_handle, fig_dir, fname_base)
% Saves .png (300 dpi), .pdf, and .fig
png_path = fullfile(fig_dir, fname_base + ".png");
pdf_path = fullfile(fig_dir, fname_base + ".pdf");
fig_path = fullfile(fig_dir, fname_base + ".fig");
try
    exportgraphics(fig_handle, png_path, 'Resolution',300);
    exportgraphics(fig_handle, pdf_path);
catch
    % fallback for older export issues
    saveas(fig_handle, png_path);
    saveas(fig_handle, pdf_path);
end
savefig(fig_handle, fig_path);
end
