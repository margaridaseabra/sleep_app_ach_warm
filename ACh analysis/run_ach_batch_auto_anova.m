function GROUP = run_ach_batch_auto_anova(sigDir, scoreDir, varargin)
% run_ach_batch_auto
% -------------------------------------------------------------------------
% Automatically:
%   - finds all ACh .mat signal files in sigDir
%   - finds all 1-Hz score CSVs in scoreDir
%   - matches them by mouse ID and condition
%   - runs ach_analysis() for each pair
%   - collects per-session metrics into GROUP and ACh_group_metrics.csv
%
% Expected file name patterns:
%   Signals (.mat), e.g.:
%       20251026-baseline-mouse12_notched50Hz_bw3.mat
%       20251001_baseline_mouse1_APP.mat
%       20251001_baseline_mouse1_APP_REM_QC.mat
%
%   Scores (.csv), e.g.:
%       20251012-baseline-mouse12-WT_scored_scores_1Hz.csv
%       20251002_baseline_mouse2_WT_scored_scores_1Hz.csv
%
%   where:
%       cond  = baseline | ambtemp | drugs (etc., no dashes inside)
%       GENO  = WT | APPPS1 (or whatever label you use)
%
% New optional arguments (name/value):
%   'showMouseIDs' (false/true)  -> plot mouse IDs next to points in group plots
%   'doANOVA'     (false/true)  -> run 2-way ANOVA (Genotype x Condition)
%
% Usage:
%   GROUP = run_ach_batch_auto('/path/to/mat_folder', ...
%                              '/path/to/csv_folder');
%
%   GROUP = run_ach_batch_auto('/path/to/mat_folder', ...
%                              '/path/to/csv_folder', ...
%                              'showMouseIDs', true, ...
%                              'doANOVA', true);
%
% Requires: ach_analysis.m on path.

    % ---------------------------------------------------------------------
    % Parse optional arguments
    % ---------------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'showMouseIDs', false, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'doANOVA',     false, @(x)islogical(x) || isnumeric(x));
    parse(p, varargin{:});
    opts = p.Results;
    opts.showMouseIDs = logical(opts.showMouseIDs);
    opts.doANOVA     = logical(opts.doANOVA);

    % ---------------------------------------------------------------------
    % 0) Directories
    % ---------------------------------------------------------------------
    if nargin < 1 || isempty(sigDir)
        sigDir = uigetdir(pwd,'Select folder with ACh .mat files');
        if sigDir == 0, error('No signal folder selected.'); end
    end
    if nargin < 2 || isempty(scoreDir)
        scoreDir = uigetdir(pwd,'Select folder with 1-Hz score CSVs');
        if scoreDir == 0, error('No score folder selected.'); end
    end

    fprintf('\nSignals dir: %s\nScores dir : %s\n', sigDir, scoreDir);

    % ---------------------------------------------------------------------
    % 1) Find and parse all score files
    % ---------------------------------------------------------------------
    scoreFiles = dir(fullfile(scoreDir, '*_scores_1Hz.csv'));
    scoreMeta  = struct([]);     % let first entry define the fields

    for k = 1:numel(scoreFiles)
        fn = scoreFiles(k).name;
        info = parse_score_filename(fn);
        if isempty(info)
            fprintf('  [WARN] Could not parse score filename: %s\n', fn);
            continue;
        end
        info.file = fullfile(scoreDir, fn);

        if isempty(scoreMeta)
            scoreMeta = info;
        else
            scoreMeta(end+1) = info; 
        end
    end

    if isempty(scoreMeta)
        error('No valid score CSVs found in %s', scoreDir);
    end

    % ---------------------------------------------------------------------
    % 2) Find and parse all signal (.mat) files, match to scores
    % ---------------------------------------------------------------------
    sigFiles = dir(fullfile(sigDir, '*.mat'));
    SESS = struct([]);           % again, let first entry define fields

    for k = 1:numel(sigFiles)
        fn = sigFiles(k).name;
        sInfo = parse_signal_filename(fn);
        if isempty(sInfo)
            fprintf('  [WARN] Could not parse signal filename: %s\n', fn);
            continue;
        end

        % Build key for matching: lower(cond) & mouse ID
        condKey  = lower(sInfo.cond);
        mouseKey = sInfo.mouse;           % number as string

        % find first score file with same cond + mouse
        idx = find(strcmpi({scoreMeta.cond}, condKey) & ...
                   strcmp({scoreMeta.mouse}, mouseKey), 1, 'first');

        if isempty(idx)
            fprintf('  [WARN] No score CSV for %s (mouse %s, cond %s)\n', ...
                    fn, mouseKey, condKey);
            continue;
        end

        sm = scoreMeta(idx);

        sess.mouse = sprintf('mouse%s', mouseKey); % e.g. 'mouse12'
        sess.geno  = sm.geno;                      % WT / APPPS1
        sess.cond  = condKey;                      % baseline / ambtemp / drugs
        sess.mat   = fullfile(sigDir, fn);
        sess.csv   = sm.file;

        if isempty(SESS)
            SESS = sess;
        else
            SESS(end+1) = sess; %#ok<AGROW>
        end
    end

    if isempty(SESS)
        error('No matched signal+score pairs found.');
    end

    fprintf('\nMatched %d sessions (signal + scores).\n', numel(SESS));

    % ---------------------------------------------------------------------
    % 3) Run ach_analysis for each session and collect metrics
    % ---------------------------------------------------------------------
    CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);   % adjust if needed

    nSess = numel(SESS);
    GROUP.sessions = SESS;
    GROUP.out      = cell(nSess,1);
    GROUP.metrics  = struct([]);

    for k = 1:nSess
        s = SESS(k);
        fprintf('\n=== Processing %s | %s | %s ===\n', ...
                s.mouse, s.geno, s.cond);

        out_prefix = sprintf('%s_%s_%s', s.mouse, s.geno, s.cond);

        OUT = ach_analysis( ...
            s.mat, s.csv, ...
            'codes', CODES, ...
            'mouse_id', s.mouse, ...
            'session', s.cond, ...
            'out_prefix', out_prefix, ...
            'verbose', true);   % set false if you prefer silence

        GROUP.out{k} = OUT;

        % ---------- Light per-session summary for group stats ----------
        GROUP.metrics(k).mouse = s.mouse;
        GROUP.metrics(k).geno  = s.geno;
        GROUP.metrics(k).cond  = s.cond;

        % NREM PSD metrics (if present)
        if isfield(OUT,'psd') && isfield(OUT.psd,'NREM') && ~isempty(OUT.psd.NREM.f)
            GROUP.metrics(k).NREM_power   = OUT.psd.NREM.band_power;
            GROUP.metrics(k).NREM_peakHz  = OUT.psd.NREM.peak_freq;
            GROUP.metrics(k).NREM_peakAmp = OUT.psd.NREM.peak_amp;
        else
            GROUP.metrics(k).NREM_power   = NaN;
            GROUP.metrics(k).NREM_peakHz  = NaN;
            GROUP.metrics(k).NREM_peakAmp = NaN;
        end

        % State-wise ACh slopes (if you implemented OUT.state_slopes)
        if isfield(OUT,'state_slopes')
            GROUP.metrics(k).slope_Wake = get_state_slope(OUT,'Wake');
            GROUP.metrics(k).slope_NREM = get_state_slope(OUT,'NREM');
            GROUP.metrics(k).slope_REM  = get_state_slope(OUT,'REM');
        else
            GROUP.metrics(k).slope_Wake = NaN;
            GROUP.metrics(k).slope_NREM = NaN;
            GROUP.metrics(k).slope_REM  = NaN;
        end

        % Wake-onset transition peaks / slopes (if present)
        idxW = [];
        if isfield(OUT,'transitions')
            idxW = find(strcmp({OUT.transitions.name},'Wake_onset'));
        end
        if ~isempty(idxW)
            Tw = OUT.transitions(idxW);
            GROUP.metrics(k).WakeOn_peak_mean  = mean(Tw.peaks,  'omitnan');
            if isfield(Tw,'slopes') && ~isempty(Tw.slopes)
                GROUP.metrics(k).WakeOn_slope_mean = mean(Tw.slopes,'omitnan');
            else
                GROUP.metrics(k).WakeOn_slope_mean = NaN;
            end
        else
            GROUP.metrics(k).WakeOn_peak_mean  = NaN;
            GROUP.metrics(k).WakeOn_slope_mean = NaN;
        end
    end

    % ---------------------------------------------------------------------
    % 4) Save group metrics table to CSV + .mat
    % ---------------------------------------------------------------------
    GROUP.metrics_tbl = struct2table(GROUP.metrics);
    csv_out = fullfile(sigDir, 'ACh_group_metrics.csv');
    mat_out = fullfile(sigDir, 'ACh_GROUP_all.mat');

    writetable(GROUP.metrics_tbl, csv_out);
    save(mat_out, 'GROUP');

    fprintf('\nSaved:\n  %s\n  %s\n', csv_out, mat_out);

    % ---------------------------------------------------------------------
    % 5) Optional: ANOVA + group plots (with optional mouse IDs)
    % ---------------------------------------------------------------------
    if opts.doANOVA || opts.showMouseIDs
        GROUP = ach_group_summary(GROUP, opts);
    end

    fprintf('Done.\n\n');
end

% =======================================================================
% --------------------- Helper functions --------------------------------
% =======================================================================

function info = parse_signal_filename(fname)
% Parse signal MAT names in either of these forms:
%   1) YYYYMMDD-cond-mouseNN_something.mat
%   2) YYYYMMDD_cond_mouseNN_something.mat
%
% Examples:
%   20251001-baseline-mouse1_notched50Hz_bw3.mat
%   20251001_baseline_mouse1_APP.mat
%   20251001_baseline_mouse1_APP_REM_QC.mat
%   20251022_drugs_mouse2_WT.mat

    % Allow both '-' and '_' as separators
    pat = '^(?<date>\d{8})[-_](?<cond>[^-_]+)[-_]mouse(?<mouse>\d+).*\.mat$';

    m = regexp(fname, pat, 'names');
    if isempty(m)
        info = [];
    else
        info = struct();
        info.date  = m.date;
        info.cond  = lower(m.cond);   % use lowercase for matching
        info.mouse = m.mouse;         % just the number as string
    end
end

function info = parse_score_filename(fname)
% Parse score CSV names in either of these forms:
%   1) YYYYMMDD-cond-mouseNN-GENO_scored_scores_1Hz.csv
%   2) YYYYMMDD_cond_mouseNN_GENO_scored_scores_1Hz.csv
%
% Examples:
%   20251012-baseline-mouse12-WT_scored_scores_1Hz.csv
%   20251002_baseline_mouse2_WT_scored_scores_1Hz.csv

    % Allow both '-' and '_' as separators between fields
    pat = ['^(?<date>\d{8})[-_](?<cond>[^-_]+)[-_]mouse(?<mouse>\d+)[-_]' ...
           '(?<geno>[^-_]+)_scored_scores_1Hz\.csv$'];

    m = regexp(fname, pat, 'names');
    if isempty(m)
        info = [];
    else
        info = struct();
        info.date  = m.date;
        info.cond  = lower(m.cond);   % lowercase for matching
        info.mouse = m.mouse;         % number as string
        info.geno  = m.geno;          % keep genotype label as-is
    end
end

function val = get_state_slope(OUT, stateName)
% Gracefully get mean slope for a state from OUT.state_slopes
    if isfield(OUT,'state_slopes') && isfield(OUT.state_slopes,stateName)
        val = OUT.state_slopes.(stateName).mean;
    else
        val = NaN;
    end
end

% =======================================================================
% --------- Group stats & plotting (ANOVA + mouse ID labels) ------------
% =======================================================================

function GROUP = ach_group_summary(GROUP, opts)
% Run 2-way ANOVA (Genotype x Condition) and make group plots
% Metrics are taken from GROUP.metrics_tbl

    T = GROUP.metrics_tbl;

    % List of scalar metrics we try to analyse
    metricNames = { ...
        'NREM_power', ...
        'NREM_peakHz', ...
        'NREM_peakAmp', ...
        'slope_Wake', ...
        'slope_NREM', ...
        'slope_REM', ...
        'WakeOn_peak_mean', ...
        'WakeOn_slope_mean'};

    GROUP.ANOVA = struct();

    for iM = 1:numel(metricNames)
        mName = metricNames{iM};
        if ~ismember(mName, T.Properties.VariableNames)
            continue;
        end

        y = T.(mName);

        % skip non-numeric or non-scalar-per-row metrics
        if ~isnumeric(y)
            continue;
        end

        % y is typically an N×1 vector; but if it is N×k, skip for now
        if size(y,2) ~= 1
            continue;
        end

        valid = ~isnan(y);
        if sum(valid) < 3
            continue; % not enough data
        end

        yv     = y(valid);
        geno   = T.geno(valid);
        cond   = T.cond(valid);
        mouse  = T.mouse(valid);

        % ---------- ANOVA ----------
        if opts.doANOVA
            try
                [p, tbl, stats] = anovan(yv, {geno, cond}, ...
                    'model', 'interaction', ...
                    'varnames', {'Genotype','Condition'}, ...
                    'display', 'off');

                GROUP.ANOVA.(mName).p     = p;
                GROUP.ANOVA.(mName).table = tbl;
                GROUP.ANOVA.(mName).stats = stats;

                pInt = NaN;
                if numel(p) >= 3, pInt = p(3); end
                fprintf('\nANOVA for %s:\n', mName);
                fprintf('  p(Genotype) = %.3g\n', p(1));
                fprintf('  p(Condition)= %.3g\n', p(2));
                if ~isnan(pInt)
                    fprintf('  p(Interaction)= %.3g\n', pInt);
                end
            catch ME
                warning('ANOVA failed for %s: %s', mName, ME.message);
            end
        end

        % ---------- Plot (with optional mouse IDs) ----------
        ach_plot_metric(yv, geno, cond, mouse, mName, opts.showMouseIDs);
    end
end

function ach_plot_metric(y, geno, cond, mouse, metricName, showMouseIDs)
% Scatter plot of per-mouse values grouped by Genotype_Condition.
% If showMouseIDs = true, labels each point with mouse ID.

    if numel(y) < 2
        return;
    end

    groupStr = strcat(geno, '_', cond);   % e.g. 'WT_baseline'
    [groupNames, ~, gIdx] = unique(groupStr, 'stable');

    x = gIdx;
    % small jitter for visibility
    x = x + 0.15*(rand(size(x))-0.5);

    figure('Name', metricName);
    scatter(x, y, 40, 'filled'); hold on;

    % group means
    for g = 1:numel(groupNames)
        yg = y(gIdx == g);
        if isempty(yg), continue; end
        plot(g, mean(yg,'omitnan'), 'kd', 'MarkerSize', 9, 'MarkerFaceColor', 'k');
    end

    xlim([0.5, numel(groupNames)+0.5]);
    set(gca, 'XTick', 1:numel(groupNames), ...
             'XTickLabel', groupNames, ...
             'XTickLabelRotation', 45);
    ylabel(metricName, 'Interpreter', 'none');
    xlabel('Genotype_Condition');
    grid on;
    title(sprintf('%s per mouse', strrep(metricName, '_', '\_')));

    if showMouseIDs
        for i = 1:numel(y)
            text(x(i)+0.02, y(i), mouse{i}, ...
                'Interpreter','none', 'FontSize', 8);
        end
    end
end
