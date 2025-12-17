function OUT = run_MA_fragmentation_index(MA, ARCH, out_dir)
% run_MA_fragmentation_index
% -------------------------------------------------------------------------
% Fragmentation analysis focused on micro-arousals (MAs) during BASELINE.
%
% REQUIRED INPUTS
%   MA   : table of micro-arousal episodes with at least:
%            - mouse        : mouse ID
%            - genotype     : 'WT' / 'APP'
%            - condition    : e.g. 'baseline'
%            - ma_on_s      : MA onset time (seconds from recording start)
%            - state_before : sleep state immediately before MA
%                             (e.g. 'NREM' or 'REM')
%
%   ARCH : sleep architecture table (like rows_overall or epoch-level) with:
%            - mouse
%            - genotype
%            - condition
%            - state
%            - total_dur_s  : duration in seconds in that state
%                             (if per-epoch, we will sum it)
%
%   out_dir : folder to save figures
%
% WHAT IT DOES (BASELINE ONLY)
%
% 1) Fragmentation index (per mouse, per state):
%    - For NREM and REM separately:
%         n_MA_state       = # MAs whose state_before == that state
%         dur_state_h      = total time in that state (hours)
%         MA_per_hour_state = n_MA_state / dur_state_h
%
%    - Plots (WT vs APP) with bars = group mean, dots = individual mice.
%    - Performs unpaired t-tests on MA_per_hour_state (APP vs WT).
%
% 2) Temporal distribution of MAs across recording:
%    - Define hour_idx = floor(ma_on_s / 3600)  (0 = first hour, etc.)
%    - For each state_before (NREM, REM) and genotype:
%          count MAs per mouse × hour_idx
%    - Plot WT vs APP time courses (mean ± SEM across mice).
%
% OUTPUT struct OUT has:
%   OUT.frag_table   : per-mouse fragmentation index (NREM & REM)
%   OUT.frag_stats   : stats table (APP vs WT) per state
%   OUT.timecourse   : table of MA counts per mouse/hour/state_before
%   OUT.files.frag_NREM : PNG for NREM frag plot
%   OUT.files.frag_REM  : PNG for REM frag plot
%   OUT.files.tc_NREM   : PNG for NREM MA time-course
%   OUT.files.tc_REM    : PNG for REM MA time-course
% -------------------------------------------------------------------------

if nargin < 3 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

OUT = struct();
OUT.files = struct();

%% ---------- 0) Normalise types & BASELINE only -------------------------
MA = MA(:,:);
ARCH = ARCH(:,:);

% turn key vars into string
MA.mouse     = string(MA.mouse);
MA.genotype  = string(MA.genotype);
MA.condition = string(MA.condition);
MA.state_before = string(MA.state_before);

ARCH.mouse     = string(ARCH.mouse);
ARCH.genotype  = string(ARCH.genotype);
ARCH.condition = string(ARCH.condition);
ARCH.state     = string(ARCH.state);

% Baseline only
MA_b    = MA( lower(strtrim(MA.condition))   == "baseline", : );
ARCH_b  = ARCH( lower(strtrim(ARCH.condition)) == "baseline", : );

if isempty(MA_b) || isempty(ARCH_b)
    warning('run_MA_fragmentation_index: no baseline rows in MA or ARCH.');
    OUT.success = false;
    OUT_MAfrag = run_MA_fragmentation_index(MA, ARCH, OUT.out_dir);
    return;
end

%% ---------- 1) Get total NREM/REM durations per mouse ------------------
% ARCH_b might be epoch-level. We summarise per mouse/genotype/state.

if ~ismember("total_dur_s", ARCH_b.Properties.VariableNames)
    error('ARCH table must have a variable named total_dur_s (seconds).');
end

Gdur = groupsummary(ARCH_b, {'mouse','genotype','state'}, 'sum','total_dur_s');
Gdur.Properties.VariableNames{end} = 'dur_s';  % rename sum_total_dur_s -> dur_s

% We only care about NREM and REM here
isNREM = Gdur.state == "NREM";
isREM  = Gdur.state == "REM";

Gdur_NREM = Gdur(isNREM,:);
Gdur_REM  = Gdur(isREM,:);

% Add hours
Gdur_NREM.dur_h = Gdur_NREM.dur_s / 3600;
Gdur_REM.dur_h  = Gdur_REM.dur_s  / 3600;

%% ---------- 2) Count MAs per mouse, separated by state_before ----------
% state_before should contain "NREM" / "REM" etc.

isNREM_MA = MA_b.state_before == "NREM";
isREM_MA  = MA_b.state_before == "REM";

% ---- NREM-carrier MAs ----
[gidN, mN, gN] = findgroups(MA_b.mouse(isNREM_MA), ...
                            MA_b.genotype(isNREM_MA));
n_MA_NREM = splitapply(@numel, MA_b.ma_on_s(isNREM_MA), gidN);

TabNREM = table(mN, gN, n_MA_NREM, ...
    'VariableNames', {'mouse','genotype','n_MA_NREM'});

% Join with NREM duration
if ~isempty(Gdur_NREM)
    TabNREM = innerjoin(TabNREM, ...
        Gdur_NREM(:,{'mouse','genotype','dur_h'}), ...
        'Keys', {'mouse','genotype'});
else
    TabNREM.dur_h = NaN(height(TabNREM),1);
end

TabNREM.MA_per_h_NREM = TabNREM.n_MA_NREM ./ max(TabNREM.dur_h,eps);
TabNREM.state = repmat("NREM", height(TabNREM),1);

% ---- REM-carrier MAs ----
[gidR, mR, gR] = findgroups(MA_b.mouse(isREM_MA), ...
                            MA_b.genotype(isREM_MA));
n_MA_REM = splitapply(@numel, MA_b.ma_on_s(isREM_MA), gidR);

TabREM = table(mR, gR, n_MA_REM, ...
    'VariableNames', {'mouse','genotype','n_MA_REM'});

if ~isempty(Gdur_REM)
    TabREM = innerjoin(TabREM, ...
        Gdur_REM(:,{'mouse','genotype','dur_h'}), ...
        'Keys', {'mouse','genotype'});
else
    TabREM.dur_h = NaN(height(TabREM),1);
end

TabREM.MA_per_h_REM = TabREM.n_MA_REM ./ max(TabREM.dur_h,eps);
TabREM.state = repmat("REM", height(TabREM),1);

% unify into a single per-mouse table
allMice = union(TabNREM.mouse, TabREM.mouse);
% To keep it simple, stack as two tables with state label:
FragN = table(TabNREM.mouse, TabNREM.genotype, ...
              TabNREM.state, TabNREM.n_MA_NREM, ...
              TabNREM.dur_h, TabNREM.MA_per_h_NREM, ...
    'VariableNames', {'mouse','genotype','state','n_MA','state_dur_h','MA_per_h'});
FragR = table(TabREM.mouse, TabREM.genotype, ...
              TabREM.state, TabREM.n_MA_REM, ...
              TabREM.dur_h, TabREM.MA_per_h_REM, ...
    'VariableNames', {'mouse','genotype','state','n_MA','state_dur_h','MA_per_h'});

FRAG = [FragN; FragR];

OUT.frag_table = FRAG;

%% ---------- 3) Stats & plots: fragmentation index (NREM / REM) ---------
states_to_show = ["NREM","REM"];

frag_stats_rows = {};
for si = 1:numel(states_to_show)
    st = states_to_show(si);
    Tst = FRAG(FRAG.state == st, :);
    if isempty(Tst)
        continue;
    end

    isWT  = Tst.genotype == "WT";
    isAPP = Tst.genotype ~= "WT";   % collapse non-WT as APP

    valsWT  = Tst.MA_per_h(isWT);
    valsAPP = Tst.MA_per_h(isAPP);

    valsWT  = valsWT(~isnan(valsWT));
    valsAPP = valsAPP(~isnan(valsAPP));

    nWT = numel(valsWT);
    nAPP = numel(valsAPP);

    if nWT >= 2 && nAPP >= 2
        [~, p_t] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
        p_rs = ranksum(valsWT, valsAPP);

        m1 = mean(valsWT);  m2 = mean(valsAPP);
        s1 = std(valsWT);   s2 = std(valsAPP);
        n1 = nWT;           n2 = nAPP;
        sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
        d  = (m2 - m1) / sp;
    else
        p_t = NaN; p_rs = NaN; d = NaN;
    end

    frag_stats_rows(end+1,:) = {char(st), nWT, nAPP, ...
        mean(valsWT,'omitnan'), mean(valsAPP,'omitnan'), ...
        p_t, p_rs, d}; %#ok<AGROW>

    % ---- plot fragmentation bar for this state ----
    figure('Color','w'); hold on;
    Xpos = [1 2];
    COL_WT      = [0.6 0.6 0.6];
    COL_APP     = [0.39 0.58 0.93];
    jit = 0.12;

    barWidth = 0.5;
    bar(Xpos(1), mean(valsWT,'omitnan'),  barWidth, ...
        'FaceColor',COL_WT,'EdgeColor','none');
    bar(Xpos(2), mean(valsAPP,'omitnan'), barWidth, ...
        'FaceColor',COL_APP,'EdgeColor','none');

    % SD or SEM? For fragmentation index I'd show SEM of per-mouse values
    semWT  = std(valsWT,'omitnan')  / max(1,sqrt(nWT));
    semAPP = std(valsAPP,'omitnan') / max(1,sqrt(nAPP));

    errorbar(Xpos(1), mean(valsWT,'omitnan'),  semWT, 'k','LineStyle','none');
    errorbar(Xpos(2), mean(valsAPP,'omitnan'), semAPP,'k','LineStyle','none');

    % jittered dots
    xw = Xpos(1) + (rand(size(valsWT))-0.5)*2*jit;
    xa = Xpos(2) + (rand(size(valsAPP))-0.5)*2*jit;
    plot(xw, valsWT,  '.', 'Color',[0.3 0.3 0.3], 'MarkerSize',12);
    plot(xa, valsAPP, '.', 'Color',[0.1 0.2 0.6], 'MarkerSize',12);

    xlim([0.5 2.5]);
    set(gca,'XTick',Xpos,'XTickLabel',{'WT','APP'},'FontSize',12);
    ylabel(sprintf('MAs per hour of %s', st));
    title(sprintf('Baseline fragmentation index: %s (MAs / h of %s)', st, st));
    set(gca,'Box','off');

    ymax = max([valsWT; valsAPP],[],'omitnan');
    if isempty(ymax) || isnan(ymax), ymax = 1; end
    ylim([0, ymax*1.4]);

    % annotate p-values
    if ~isnan(p_t)
        txt = sprintf('t-test p=%.3f; ranksum p=%.3f; d=%.2f', p_t, p_rs, d);
        text(1.5, ymax*1.25, txt, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', ...
             'FontSize',10);
    end

    % save
    fname = sprintf('baseline_MA_fragmentation_%s_WTvsAPP.png', lower(st));
    fpath = fullfile(out_dir, fname);
    saveas(gcf, fpath);
    if st == "NREM"
        OUT.files.frag_NREM = fpath;
    elseif st == "REM"
        OUT.files.frag_REM  = fpath;
    end

    fprintf('✅ Fragmentation plot for %s saved to %s\n', st, fpath);
end

% stats table
if ~isempty(frag_stats_rows)
    frag_stats = cell2table(frag_stats_rows, ...
        'VariableNames', {'state','nWT','nAPP', ...
                          'mean_WT_MA_per_h','mean_APP_MA_per_h', ...
                          'p_ttest','p_ranksum','Cohen_d'});
else
    frag_stats = table();
end
OUT.frag_stats = frag_stats;

if ~isempty(frag_stats)
    fprintf('\nFragmentation index (MAs per hour of state): WT vs APP\n');
    disp(frag_stats);
end

%% ---------- 4) Temporal distribution of MAs across hours ---------------
% We now look at *when* MAs occur during the baseline recording.

MA_tc = MA_b;
% hour_idx: 0 = first hour [0,3600), 1 = [3600,7200), etc.
MA_tc.hour_idx = floor(MA_tc.ma_on_s / 3600);

% We'll keep only NREM/REM carrier states for the time-course
isNR = MA_tc.state_before == "NREM";
isRM = MA_tc.state_before == "REM";

states_tc = ["NREM","REM"];
OUT.timecourse = struct();

for si = 1:numel(states_tc)
    st = states_tc(si);
    if st == "NREM"
        Mst = MA_tc(isNR,:);
    else
        Mst = MA_tc(isRM,:);
    end
    if isempty(Mst), continue; end

    [gid_tc, m_tc, g_tc, h_tc] = findgroups(Mst.mouse, Mst.genotype, Mst.hour_idx);
    n_MA_hr = splitapply(@numel, Mst.ma_on_s, gid_tc);

    T_tc = table(m_tc, g_tc, h_tc, n_MA_hr, ...
        'VariableNames', {'mouse','genotype','hour_idx','n_MA'});

    OUT.timecourse.(lower(st)) = T_tc;

    % get list of hours in this state
    all_hours = unique(T_tc.hour_idx);
    all_hours = sort(all_hours);

    meanWT  = nan(size(all_hours));
    semWT   = nan(size(all_hours));
    meanAPP = nan(size(all_hours));
    semAPP  = nan(size(all_hours));

    for k = 1:numel(all_hours)
        hr = all_hours(k);
        maskH = T_tc.hour_idx == hr;

        valsWT  = T_tc.n_MA(maskH & T_tc.genotype == "WT");
        valsAPP = T_tc.n_MA(maskH & T_tc.genotype ~= "WT");

        valsWT  = valsWT(~isnan(valsWT));
        valsAPP = valsAPP(~isnan(valsAPP));

        if ~isempty(valsWT)
            meanWT(k) = mean(valsWT);
            semWT(k)  = std(valsWT) / max(1,sqrt(numel(valsWT)));
        end
        if ~isempty(valsAPP)
            meanAPP(k) = mean(valsAPP);
            semAPP(k)  = std(valsAPP) / max(1,sqrt(numel(valsAPP)));
        end
    end

    % ---- plot time-course ----
    figure('Color','w'); hold on;
    COL_WT  = [0.6 0.6 0.6];
    COL_APP = [0.39 0.58 0.93];

    if any(~isnan(meanWT))
        errorbar(all_hours, meanWT, semWT, '-o', ...
            'Color',COL_WT, 'MarkerFaceColor',COL_WT, ...
            'LineWidth',1.2);
    end
    if any(~isnan(meanAPP))
        errorbar(all_hours, meanAPP, semAPP, '-o', ...
            'Color',COL_APP,'MarkerFaceColor',COL_APP, ...
            'LineWidth',1.2);
    end

    xlabel('Hour from recording start');
    ylabel(sprintf('MAs per mouse (carrier=%s)', st));
    title(sprintf('Baseline temporal distribution of MAs (carrier=%s)', st));
    legend({'WT','APP'}, 'Location','northoutside', 'Orientation','horizontal');
    set(gca,'Box','off','FontSize',12);
    xlim([min(all_hours)-0.5, max(all_hours)+0.5]);

    fname_tc = sprintf('baseline_MA_timecourse_%s_WTvsAPP.png', lower(st));
    fpath_tc = fullfile(out_dir, fname_tc);
    saveas(gcf, fpath_tc);

    if st == "NREM"
        OUT.files.tc_NREM = fpath_tc;
    elseif st == "REM"
        OUT.files.tc_REM  = fpath_tc;
    end

    fprintf('✅ MA time-course plot for carrier=%s saved to %s\n', st, fpath_tc);
end

OUT.success = true;
end
