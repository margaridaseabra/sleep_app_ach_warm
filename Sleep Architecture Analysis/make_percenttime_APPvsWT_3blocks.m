function OUT = make_percenttime_APPvsWT_3blocks(PERHOUR, out_dir, ...
                                                states_to_use, conds_to_use, varargin)
% make_percenttime_APPvsWT_3blocks
% -------------------------------------------------------------------------
% From per-hour sleep architecture (OUT.per_hour from
% run_group_sleep_architecture), compute % time in each state for:
%   - 3 time blocks: 0–3 h, 3–6 h, washout (>= 6 h)
%   - selected conditions (e.g. baseline, ambtemp, drug)
%   - genotypes: WT vs APP (any non-WT collapsed to APP)
%
% For each STATE in states_to_use, makes ONE figure:
%   X-axis : 3 blocks (0–3 h, 3–6 h, washout)
%   Within each block: bars for every (condition x genotype) combo
%                      e.g. WT baseline, WT ambtemp, APP baseline, APP ambtemp
%
% Stats:
%   1) Global 3-way ANOVA (block x condition x genotype), printed top-left.
%   2) Per-block 2-way ANOVAs (condition x genotype), printed under each block:
%         "<block>: c p=.., g p=.., int p=.."
%
% Visual:
%   - 3–6 h block is shaded and labelled "manipulation" in every state plot.
%
% INPUT
%   PERHOUR      : table = OUT.per_hour from run_group_sleep_architecture
%   out_dir      : folder for saving figures
%   states_to_use: string array or cellstr, e.g. ["WK","NREM","REM"]
%   conds_to_use : string array or cellstr, e.g.
%                     ["baseline","ambtemp"]
%                     ["baseline","ambtemp","drug"]
%
% OPTIONAL ('Name',value):
%   'showIDs' : true/false, show mouse IDs next to dots (default = false)
%   'minN'    : minimum total N for ANOVA (default = 2)
%
% OUTPUT
%   OUT.success
%   OUT.state_stats(s):
%       .state
%       .p_block, .p_cond, .p_geno (global)
%       .block(b).name
%       .block(b).p_cond, .p_geno, .p_inter (per block)
% -------------------------------------------------------------------------

% ---------- defaults ----------
if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

if nargin < 3 || isempty(states_to_use)
    states_to_use = ["WK","MA","NREM","REM"];
end
states_to_use = string(states_to_use(:)).';

if nargin < 4 || isempty(conds_to_use)
    conds_to_use = ["baseline","ambtemp","drug"];
end
conds_to_use = lower(string(conds_to_use(:)).');

% ---------- optional args ----------
p = inputParser;
addParameter(p,'showIDs',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minN',2,@(x)isscalar(x)&&x>=1);
parse(p, varargin{:});
showIDs = p.Results.showIDs;
minN    = p.Results.minN;

PH = PERHOUR;

% ---------- 1) Check required columns ----------
mustHave = {'hour_idx','dur_s','state','mouse','genotype','condition'};
if ~all(ismember(mustHave, PH.Properties.VariableNames))
    error('PERHOUR table must contain columns: %s', strjoin(mustHave,', '));
end

% ---------- 2) Normalise text columns ----------
PH.state     = string(PH.state);
PH.mouse     = string(PH.mouse);
PH.genotype  = string(PH.genotype);
PH.condition = lower(string(strtrim(PH.condition)));
PH.condition(PH.condition=="drugs") = "drug";   % unify name

% ---------- 3) Create 3h blocks from hour_idx ----------
h = double(PH.hour_idx);   % 0,1,2,...
blk = strings(height(PH),1);
blk(h >= 0 & h < 3) = "0-3h";
blk(h >= 3 & h < 6) = "3-6h";
blk(h >= 6)         = "washout";
PH.block3h = blk;

% ---------- 4) Filter by conditions & states ----------
PH = PH(ismember(PH.condition, conds_to_use) & ...
        ismember(PH.state, states_to_use), :);

if isempty(PH)
    warning('No rows after filtering by condition/state.');
    OUT = struct('success',false,'msg','no data');
    return;
end

% Collapse any non-WT genotype to APP
geno = PH.genotype;
geno(geno ~= "WT") = "APP";
PH.geno_group = geno;

% Keep only canonical blocks in fixed order
blk_order = ["0-3h","3-6h","washout"];
PH.block3h = string(PH.block3h);
PH = PH(ismember(PH.block3h, blk_order), :);
if isempty(PH)
    warning('No rows within 0–3 h / 3–6 h / washout blocks.');
    OUT = struct('success',false,'msg','no block data');
    return;
end

% ---------- 5) Aggregate to % time per mouse × geno × cond × block × state ----------
% Sum duration per mouse/genotype/cond/block/state
[gid, m, g, c, b, st] = findgroups( ...
    PH.mouse, PH.geno_group, PH.condition, PH.block3h, PH.state);

state_dur = splitapply(@(x) sum(double(x),'omitnan'), PH.dur_s, gid);

Tstate = table(m, g, c, b, st, state_dur, ...
               'VariableNames', {'mouse','geno','cond','block','state','state_dur_s'});

% Total time per mouse × geno × cond × block (sum across states)
[gid2, m2, g2, c2, b2] = findgroups(Tstate.mouse, Tstate.geno, ...
                                    Tstate.cond, Tstate.block);
block_dur = splitapply(@(x) sum(double(x),'omitnan'), Tstate.state_dur_s, gid2);

Tblock = table(m2, g2, c2, b2, block_dur, ...
               'VariableNames', {'mouse','geno','cond','block','block_dur_s'});

% Join back to get % time
Tall = outerjoin(Tstate, Tblock, ...
                 'Keys', {'mouse','geno','cond','block'}, ...
                 'MergeKeys', true);

Tall.pct_time = Tall.state_dur_s ./ max(Tall.block_dur_s,1) * 100;

% ---------- 6) Set orders ----------
cond_order = conds_to_use;
geno_order = ["WT","APP"];
nStates    = numel(states_to_use);
nBlocks    = numel(blk_order);
nConds     = numel(cond_order);
nGen       = numel(geno_order);

% colours: base grey for WT, blue for APP, varied per condition
base_WT  = [0.6 0.6 0.6];
base_APP = [0.39 0.58 0.93];
scale_f  = linspace(1.1, 0.7, nConds);  % brightness variation

COL_BAR = zeros(nConds, 2, 3); % (cond, geno, rgb)
for j = 1:nConds
    COL_BAR(j,1,:) = min(base_WT  * scale_f(j), 1);  % WT
    COL_BAR(j,2,:) = min(base_APP * scale_f(j), 1);  % APP
end

OUT = struct();
OUT.success      = true;
OUT.state_stats  = struct([]);

% ================== 7) Loop over states ==========================
for sIdx = 1:nStates
    st_name = states_to_use(sIdx);
    Ts = Tall(Tall.state == st_name, :);
    if isempty(Ts), continue; end

    fig = figure('Color','w', ...
                 'Name',sprintf('%s – percent time',st_name));
    ax  = axes(fig); hold(ax,'on');

    % Collect data for ANOVA
    y_all        = [];
    blk_idx_all  = [];
    cond_idx_all = [];
    geno_idx_all = [];

    % X positions: 3 blocks; within each block, bars for cond×geno
    nBarsPerBlock = nConds * nGen;
    groupWidth    = 0.8;
    barWidth      = groupWidth / nBarsPerBlock;

    legend_entries = {};
    legend_handles = [];

    for bIdx = 1:nBlocks
        blk = blk_order(bIdx);

        for cIdx = 1:nConds
            cond = cond_order(cIdx);

            for gIdx = 1:nGen
                geno = geno_order(gIdx);

                % bar index within block
                bNum = (cIdx-1)*nGen + gIdx;
                x_center = bIdx + (bNum - (nBarsPerBlock+1)/2) * barWidth;

                mask = (Ts.block == blk) & (Ts.cond == cond) & (Ts.geno == geno);
                vals = Ts.pct_time(mask);
                vals = vals(~isnan(vals));

                if isempty(vals)
                    continue;
                end

                mu = mean(vals);
                se = std(vals) / max(1,sqrt(numel(vals)));

                colBar = squeeze(COL_BAR(cIdx,gIdx,:)).';
                hBar = bar(ax, x_center, mu, barWidth, ...
                           'FaceColor', colBar, 'EdgeColor','none');

                errorbar(ax, x_center, mu, se, 'k', ...
                         'LineStyle','none','LineWidth',1);

                % scatter individual mice
                x_jit = x_center + (rand(size(vals))-0.5)*barWidth*0.6;
                plot(ax, x_jit, vals, '.', ...
                     'Color',[0.2 0.2 0.2],'MarkerSize',10);

                if showIDs
                    ids = Ts.mouse(mask);
                    for q = 1:numel(vals)
                        text(ax, x_jit(q), vals(q), char(ids(q)), ...
                             'Rotation',45, ...
                             'HorizontalAlignment','left', ...
                             'VerticalAlignment','bottom', ...
                             'FontSize',7, ...
                             'Color',[0.2 0.2 0.2]);
                    end
                end

                % build legend entries for this figure
                legLabel = sprintf('%s %s', geno, cond);
                legend_entries{end+1} = legLabel; %#ok<AGROW>
                legend_handles(end+1) = hBar;    %#ok<AGROW>

                % collect for ANOVA
                nVals = numel(vals);
                y_all        = [y_all; vals];                %#ok<AGROW>
                blk_idx_all  = [blk_idx_all;  repmat(bIdx, nVals, 1)]; %#ok<AGROW>
                cond_idx_all = [cond_idx_all; repmat(cIdx, nVals, 1)]; %#ok<AGROW>
                geno_idx_all = [geno_idx_all; repmat(gIdx, nVals, 1)]; %#ok<AGROW>
            end
        end
    end

    % Axes cosmetics before shading
    xlim(ax,[0.5, nBlocks+0.5]);
    set(ax, 'XTick', 1:nBlocks, 'XTickLabel', blk_order);
    ylabel(ax,'% time in state');
    xlabel(ax,'Time block');
    title(ax, sprintf('%% time in %s – APP vs WT (%s)', ...
                      st_name, strjoin(conds_to_use,', ')));
    box(ax,'off');

    % --- Shade the 3–6 h "manipulation" block (block index 2) ---
    yl   = ylim(ax);
    ymin = yl(1);
    ymax = yl(2);
    bManip = 2;  % "3-6h"
    if nBlocks >= 2
        x1 = bManip - 0.5;
        x2 = bManip + 0.5;
        patch(ax, [x1 x2 x2 x1], [ymin ymin ymax ymax], ...
              [0.9 0.9 1.0], ...     % light bluish shade
              'FaceAlpha', 0.15, ...
              'EdgeColor','none', ...
              'HandleVisibility','off');
        % label "manipulation" – appears on every state's figure
        text(ax, bManip, ymax - 0.05*(ymax-ymin), 'manipulation', ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','top', ...
             'FontSize',9, ...
             'Color',[0.2 0.2 0.5], ...
             'FontWeight','bold');
    end

    % Legend (for every state figure)
    if ~isempty(legend_handles)
        [~, ia] = unique(legend_entries, 'stable');
        legend(ax, legend_handles(ia), legend_entries(ia), ...
               'Location','northoutside', ...
               'Orientation','horizontal', ...
               'Box','off');
    end

    % ----------------- Global 3-way ANOVA for this state ----------------
    ok = ~isnan(y_all);
    yv = y_all(ok);
    bf = blk_idx_all(ok);
    cf = cond_idx_all(ok);
    gf = geno_idx_all(ok);

    OUT.state_stats(sIdx).state = st_name;

    if numel(yv) >= minN
        [p, tbl] = anovan(yv, {bf, cf, gf}, ...
                          'model','full', ...
                          'display','off', ...
                          'varnames', {'block','cond','geno'});
        p_block = p(1);
        p_cond  = p(2);
        p_geno  = p(3);

        txt = sprintf(['Global ANOVA: block p=%.3g  |  cond p=%.3g  |  geno p=%.3g\n' ...
                       'n = %d values'], ...
                      p_block, p_cond, p_geno, numel(yv));
        text(ax, 0.01, 0.98, txt, ...
             'Units','normalized', ...
             'HorizontalAlignment','left', ...
             'VerticalAlignment','top', ...
             'FontSize',8);

        OUT.state_stats(sIdx).p_block     = p_block;
        OUT.state_stats(sIdx).p_cond      = p_cond;
        OUT.state_stats(sIdx).p_geno      = p_geno;
        OUT.state_stats(sIdx).anova_table = tbl;
    else
        OUT.state_stats(sIdx).p_block = NaN;
        OUT.state_stats(sIdx).p_cond  = NaN;
        OUT.state_stats(sIdx).p_geno  = NaN;
        OUT.state_stats(sIdx).anova_table = [];
        text(ax, 0.01, 0.98, 'Not enough data for global ANOVA', ...
             'Units','normalized', ...
             'HorizontalAlignment','left', ...
             'VerticalAlignment','top', ...
             'FontSize',8);
    end

    % ------------ Block-wise ANOVAs (cond x geno within each block) -----
    OUT.state_stats(sIdx).block = struct([]);
    p_block_cond   = NaN(1,nBlocks);
    p_block_geno   = NaN(1,nBlocks);
    p_block_inter  = NaN(1,nBlocks);

    for bIdx = 1:nBlocks
        maskB = ok & (bf == bIdx);
        yB    = y_all(maskB);
        cB    = cond_idx_all(maskB);
        gB    = geno_idx_all(maskB);

        OUT.state_stats(sIdx).block(bIdx).name = blk_order(bIdx);

        if numel(yB) < minN || numel(unique(cB)) < 2 || numel(unique(gB)) < 2
            OUT.state_stats(sIdx).block(bIdx).p_cond  = NaN;
            OUT.state_stats(sIdx).block(bIdx).p_geno  = NaN;
            OUT.state_stats(sIdx).block(bIdx).p_inter = NaN;
            continue;
        end

        % full 2-way ANOVA: cond x geno (with interaction)
        [pB, tblB] = anovan(yB, {cB, gB}, ...
                            'model','full', ...
                            'display','off', ...
                            'varnames', {'cond','geno'});

        p_block_cond(bIdx)  = pB(1);
        p_block_geno(bIdx)  = pB(2);
        p_block_inter(bIdx) = pB(3);

        OUT.state_stats(sIdx).block(bIdx).p_cond      = pB(1);
        OUT.state_stats(sIdx).block(bIdx).p_geno      = pB(2);
        OUT.state_stats(sIdx).block(bIdx).p_inter     = pB(3);
        OUT.state_stats(sIdx).block(bIdx).anova_table = tblB;
    end

    % -------- Show block-wise stats for ALL blocks on the figure --------
    yl2   = ylim(ax);
    ymin2 = yl2(1);
    ymax2 = yl2(2);
    for bIdx = 1:nBlocks
        if isnan(p_block_cond(bIdx)) || isnan(p_block_geno(bIdx))
            continue;
        end
        x0 = bIdx;
        txtB = sprintf('%s: c p=%.2g, g p=%.2g, int p=%.2g', ...
                       blk_order(bIdx), ...
                       p_block_cond(bIdx), ...
                       p_block_geno(bIdx), ...
                       p_block_inter(bIdx));
        col = [0.2 0.2 0.2];
        if bIdx == 2   % highlight manipulation block p values
            col = [0.2 0.2 0.5];
        end
        text(ax, x0, ymin2 + 0.03*(ymax2-ymin2), txtB, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', ...
             'FontSize',7, ...
             'Color', col);
    end

    % Save figure
    fname = sprintf('PercentTime_%s_%s.png', ...
                    char(st_name), strjoin(conds_to_use,'_'));
    fname = strrep(fname,'/','-');  % safety
    saveas(fig, fullfile(out_dir, fname));
end
end
