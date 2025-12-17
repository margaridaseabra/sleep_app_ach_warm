function make_state_percent_by_cond_genotype(PERHOUR, out_dir, varargin)
% make_state_percent_by_cond_genotype
% -------------------------------------------------------------------------
% For each state (WK, NREM, REM), plot:
%   "% time in state (of WK+NREM+REM)" per CONDITION,
%   with 2 bars per condition: WT vs APP.
%
% Layout:
%   [ WT  APP ]   <gap>   [ WT  APP ]   <gap>   [ WT  APP ]
%    baseline              ambient              drug
%
% SPECIAL FOR THIS PLOT:
%   For "ambient warming" and "drug" conditions we EXCLUDE
%   the first 1.5 h (5400 s) of recording (based on hour_start_s).
%
% Uses PERHOUR table created by run_group_sleep_architecture:
%   hour_idx, hour_start_s, dur_s, state, condition, mouse, genotype, ...
%
% OPTIONS (name-value):
%   'cond_order' : cellstr desired order of conditions
%                  (default tries {'baseline','ambtemp','drugs'})
%   'crop_first_sec' : scalar, number of seconds to drop for
%                      non-baseline conditions (default 5400 = 1.5 h)
%   'COL_WT' : [1x3] RGB for WT bars (default [0.6 0.6 0.6])
%   'COL_APP': [1x3] RGB for APP bars (default [0.2 0.45 0.9])

    if nargin < 2 || isempty(out_dir)
        out_dir = pwd;
    end
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    ip = inputParser;
    addParameter(ip,'cond_order', {'baseline','ambtemp','drugs'}, @(c)iscellstr(c));
    addParameter(ip,'crop_first_sec', 5400, @(x)isscalar(x)&&x>=0);
    addParameter(ip,'COL_WT',  [0.6 0.6 0.6], @(x)isnumeric(x)&&numel(x)==3);
    addParameter(ip,'COL_APP', [0.2 0.45 0.9], @(x)isnumeric(x)&&numel(x)==3);
    parse(ip, varargin{:});
    OPT = ip.Results;

    % ---------- basic guard ----------
    if isempty(PERHOUR) || height(PERHOUR)==0
        warning('make_state_percent_by_cond_genotype: PERHOUR is empty, skipping plot.');
        return;
    end

    % Work on a copy so we can crop without side-effects
    T = PERHOUR;

    % ---------- define which conditions get cropping ----------
    all_conds   = unique(T.condition); %#ok<NASGU>
    is_baseline = strcmpi(T.condition, 'baseline');
    is_manip    = ~is_baseline;   % ambient warming, drug, etc.

    % Crop for manipulation conditions: remove rows with hour_start_s < 1.5 h
    mask_keep = true(height(T),1);
    mask_keep(is_manip & T.hour_start_s < OPT.crop_first_sec) = false;
    T = T(mask_keep, :);

    % ---------- restrict to main sleep/wake states ----------
    main_states = {'WK','NREM','REM'};
    T = T(ismember(T.state, main_states), :);
    if isempty(T)
        warning('No WK/NREM/REM rows left after cropping – nothing to plot.');
        return;
    end

    % ---------- aggregate total duration per mouse x genotype x condition x state ----------
    G = groupsummary(T, {'mouse','genotype','condition','state'}, 'sum','dur_s');
    % G has variables: mouse, genotype, condition, state, GroupCount, sum_dur_s

    % Pivot to wide format: one row per mouse/genotype/condition, columns = states
    W = unstack(G, 'sum_dur_s', 'state');
    % After this, W has columns:
    %   mouse, genotype, condition, WK, NREM, REM (and maybe MA if present)

    % Ensure absent state columns exist as zeros
    for s = 1:numel(main_states)
        st = main_states{s};
        if ~ismember(st, W.Properties.VariableNames)
            W.(st) = 0;
        end
    end

    % ---------- compute percentages per row ----------
    total_main = W.WK + W.NREM + W.REM;
    total_main(total_main == 0) = NaN;  % avoid division-by-zero

    W.pct_WK   = 100 * W.WK   ./ total_main;
    W.pct_NREM = 100 * W.NREM ./ total_main;
    W.pct_REM  = 100 * W.REM  ./ total_main;

    % ---------- define condition order ----------
    cond_order = {};
    for i = 1:numel(OPT.cond_order)
        if any(strcmpi(W.condition, OPT.cond_order{i}))
            cond_order{end+1} = OPT.cond_order{i}; %#ok<AGROW>
        end
    end
    % add any other conditions that weren't in cond_order
    others = setdiff(cellstr(unique(W.condition)), cond_order, 'stable');
    cond_order = [cond_order, others];

    if isempty(cond_order)
        warning('No conditions found in PERHOUR; skipping plot.');
        return;
    end

    genotypes = {'WT','APP'};
    nCond = numel(cond_order);
    nGen  = numel(genotypes);

    % Precompute x positions: [1 2]  [4 5]  [7 8] ...
    x_bar = zeros(nCond, nGen);
    for c = 1:nCond
        base = (c-1)*(nGen+1);
        x_bar(c,:) = base + (1:nGen);
    end
    xticks_vals = mean(x_bar, 2);

    % ---------- helper function for a single state ----------
    function make_one_state_plot(state_label, pct_field)
        fh = figure('Color','w','Name', ...
            sprintf('State %% %s by condition/genotype',state_label));
        hold on;

        means = nan(nCond, nGen);
        sems  = nan(nCond, nGen);

        % bars + error bars + per-mouse dots with labels
        for c = 1:nCond
            cond_name = cond_order{c};
            for g = 1:nGen
                geno = genotypes{g};
                rows = strcmpi(W.condition, cond_name) & strcmpi(W.genotype, geno);

                vals   = W.(pct_field)(rows);
                miceID = W.mouse(rows);

                vals   = vals(~isnan(vals));
                miceID = miceID(~isnan(vals));  % keep same mask (NaNs removed)

                if isempty(vals)
                    continue;
                end

                means(c,g) = mean(vals);
                sems(c,g)  = std(vals) / sqrt(numel(vals));

                % choose colour
                switch upper(geno)
                    case 'WT'
                        col = OPT.COL_WT;
                    otherwise
                        col = OPT.COL_APP;
                end

                xb = x_bar(c,g);

                % bar + errorbar
                bar(xb, means(c,g), 'FaceColor', col, 'EdgeColor','k');
                if ~isnan(sems(c,g))
                    errorbar(xb, means(c,g), sems(c,g), ...
                        'k','LineStyle','none','LineWidth',1);
                end

                % ---- overlay per-mouse dots + mouse IDs ----
                jitter = (rand(size(vals))-0.5)*0.20;
                x_pts  = xb + jitter;

                scatter(x_pts, vals, 28, col, 'filled', ...
                        'MarkerFaceAlpha',0.8, ...
                        'MarkerEdgeColor','k', 'MarkerEdgeAlpha',0.6);

                % mouse ID labels
                for k = 1:numel(vals)
                    text(x_pts(k)+0.05, vals(k), char(miceID(k)), ...
                         'FontSize',8, ...
                         'Color', col, ...
                         'HorizontalAlignment','left', ...
                         'VerticalAlignment','middle');
                end
            end
        end

        % cosmetics
        xlim([0, max(x_bar(:)) + 1]);
        ylim([0 100]);
        set(gca,'XTick', xticks_vals, 'XTickLabel', cond_order, 'FontSize',10);
        ylabel(sprintf('Time in %s (%% of WK+NREM+REM)', state_label));
        xlabel('Condition');
        title(sprintf('%% Time in %s by condition and genotype', state_label), ...
              'Interpreter','none');

        % legend patches
        p1 = patch(NaN,NaN,OPT.COL_WT,'EdgeColor','k');
        p2 = patch(NaN,NaN,OPT.COL_APP,'EdgeColor','k');
        legend([p1 p2], {'WT','APP'}, 'Location','best','Box','off');

        box on; grid on;

        % save
        fname = fullfile(out_dir, ...
            sprintf('state_percent_%s_by_condition_genotype.png', ...
                    lower(state_label)));
        exportgraphics(fh, fname, 'Resolution',300);
        fprintf('Saved %s\n', fname);
    end

    % ---------- make one figure per state ----------
    make_one_state_plot('WK',   'pct_WK');
    make_one_state_plot('NREM','pct_NREM');
    make_one_state_plot('REM', 'pct_REM');
end
