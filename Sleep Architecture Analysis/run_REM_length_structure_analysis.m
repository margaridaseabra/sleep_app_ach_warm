function OUT = run_REM_length_structure_analysis(input_dir, varargin)
% run_REM_length_structure_analysis
% -------------------------------------------------------------------------
% 1. For ALL recordings in input_dir:
%       - run sleep_architecture_from_scores to get Runs table
%       - extract all REM bouts (start_s, end_s, dur_s)
%       - attach file/date/condition/mouse/genotype (parsed from filename)
%
% 2. Build GLOBAL REM duration distribution (all conditions) and define:
%       short / medium / long REM using quantiles (default [33 66]).
%
% 3. BASELINE-ONLY analysis:
%       - per mouse: count short / medium / long REM bouts
%       - summarize APP vs WT
%       - make bar plot: #REM bouts per type (Short/Medium/Long), WT vs APP
%
% 4. REM CLUSTER STRUCTURE (baseline only):
%       - clusters of REM bouts separated by gaps <= cluster_gap_sec
%       - cluster types: ShortOnly, ShortLong, LongOnly, MediumOnly
%       - per long REM: n_short_before_long
%       - per mouse summary of cluster structure.
% -------------------------------------------------------------------------

p = inputParser;
addRequired(p,'input_dir',@ischar);
addParameter(p,'pattern','*_scores_1Hz.csv',@ischar);
addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
addParameter(p,'quantiles',[33 66],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'cluster_gap_sec',600,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'out_dir','',@ischar);
parse(p, input_dir, varargin{:});
S = p.Results;
C = S.codes;

assert(isfolder(S.input_dir), 'Input folder not found: %s', S.input_dir);

if isempty(S.out_dir)
    S.out_dir = fullfile(S.input_dir, 'rem_analysis');
end
if ~isfolder(S.out_dir)
    mkdir(S.out_dir);
end

if exist('sleep_architecture_from_scores','file') ~= 2
    error(['sleep_architecture_from_scores.m is not on the MATLAB path.\n' ...
           'Add it, e.g.: addpath(genpath(''%s''))'], pwd);
end

% =======================
% 1) DISCOVER FILES
% =======================
F = dir(fullfile(S.input_dir, '**', S.pattern));
F = F(~[F.isdir]);
assert(~isempty(F), 'No files matched "%s" under %s', S.pattern, S.input_dir);

ALL_REM = table();

% =======================
% 2) LOOP FILES, GET REM
% =======================
for i = 1:numel(F)
    csv_path = fullfile(F(i).folder, F(i).name);
    info = parse_info_from_filename_relaxed(F(i).name);  % <--- NEW PARSER

    if ~info.ok
        warning('Filename not recognized pattern, skipping: %s', F(i).name);
        continue;
    end

    try
        OUT_i = sleep_architecture_from_scores(csv_path, ...
            'codes', C, ...
            'includeMA', true, ...
            'ma_thresh_sec', 15, ...
            'reclassify_short_wake_to_MA', true, ...
            'out_prefix','', ...
            'out_dir','', ...
            'write_log', false, ...
            'verbose', false);
    catch ME
        warning('Failed file %s: %s', F(i).name, ME.message);
        continue;
    end

    Runs = OUT_i.runs;
    if isempty(Runs); continue; end

    rem_mask = (Runs.state_code == C.REM);
    if ~any(rem_mask); continue; end

    R = Runs(rem_mask, :);
    n = height(R);

    R.file      = repmat(string(F(i).name),      n,1);
    R.date      = repmat(string(info.date),      n,1);
    R.condition = repmat(string(info.condition), n,1);
    R.mouse     = repmat(string(info.mouse),     n,1);
    R.genotype  = repmat(string(info.genotype),  n,1);

    ALL_REM = [ALL_REM; R]; %#ok<AGROW>
end

if isempty(ALL_REM)
    error('No REM bouts found in any file (check codes and filename parsing).');
end

% ensure dur_s
if ~ismember('dur_s', ALL_REM.Properties.VariableNames)
    ALL_REM.dur_s = ALL_REM.end_s - ALL_REM.start_s + 1;
end

% ================================
% 3) GLOBAL SHORT/MEDIUM/LONG DEF
% ================================
rem_durs = double(ALL_REM.dur_s);
rem_durs = rem_durs(~isnan(rem_durs) & rem_durs > 0);

q = prctile(rem_durs, S.quantiles);
q1 = q(1); q2 = q(2);

REM_THRESH = struct();
REM_THRESH.quantiles = S.quantiles;
REM_THRESH.q1_sec    = q1;
REM_THRESH.q2_sec    = q2;

rem_type = strings(height(ALL_REM),1);
rem_type(ALL_REM.dur_s <  q1) = "Short";
rem_type(ALL_REM.dur_s >= q1 & ALL_REM.dur_s < q2) = "Medium";
rem_type(ALL_REM.dur_s >= q2) = "Long";
ALL_REM.rem_type = rem_type;

all_rem_csv = fullfile(S.out_dir, 'ALL_REM_with_types.csv');
writetable(ALL_REM, all_rem_csv);

% ================================
% 4) BASELINE-ONLY COUNTS (APP vs WT)
% ================================
is_baseline = lower(strtrim(ALL_REM.condition)) == "baseline";
BASE_REM = ALL_REM(is_baseline, :);

COUNTS_PER_MOUSE = table();
COUNTS_BY_GENO   = table();
counts_csv = ''; counts_geno_csv = ''; rem_counts_png = '';

if ~isempty(BASE_REM)
    mice = unique(BASE_REM.mouse);
    for m = 1:numel(mice)
        mID = mice(m);
        maskM = BASE_REM.mouse == mID;
        Gm = BASE_REM(maskM,:);

        geno = unique(Gm.genotype);
        geno = geno(1);

        nShort  = nnz(Gm.rem_type == "Short");
        nMedium = nnz(Gm.rem_type == "Medium");
        nLong   = nnz(Gm.rem_type == "Long");

        COUNTS_PER_MOUSE = [COUNTS_PER_MOUSE; ...
            table(mID, geno, nShort, nMedium, nLong, ...
                  'VariableNames', {'mouse','genotype','nShort','nMedium','nLong'})]; %#ok<AGROW>
    end

    counts_csv = fullfile(S.out_dir, 'baseline_REM_counts_per_mouse.csv');
    writetable(COUNTS_PER_MOUSE, counts_csv);

    genos = unique(COUNTS_PER_MOUSE.genotype);
    for g = 1:numel(genos)
        gg = genos(g);
        maskG = COUNTS_PER_MOUSE.genotype == gg;
        n_mice = nnz(maskG);

        meanShort  = mean(COUNTS_PER_MOUSE.nShort(maskG));
        meanMedium = mean(COUNTS_PER_MOUSE.nMedium(maskG));
        meanLong   = mean(COUNTS_PER_MOUSE.nLong(maskG));

        semShort  = std(COUNTS_PER_MOUSE.nShort(maskG))  / sqrt(n_mice);
        semMedium = std(COUNTS_PER_MOUSE.nMedium(maskG)) / sqrt(n_mice);
        semLong   = std(COUNTS_PER_MOUSE.nLong(maskG))   / sqrt(n_mice);

        COUNTS_BY_GENO = [COUNTS_BY_GENO; ...
            table(gg, n_mice, meanShort, semShort, ...
                  meanMedium, semMedium, ...
                  meanLong, semLong, ...
                  'VariableNames', {'genotype','n_mice', ...
                                    'meanShort','semShort', ...
                                    'meanMedium','semMedium', ...
                                    'meanLong','semLong'})]; %#ok<AGROW>
    end

    counts_geno_csv = fullfile(S.out_dir, 'baseline_REM_counts_by_genotype.csv');
    writetable(COUNTS_BY_GENO, counts_geno_csv);

    % ---- bar plot of #short/medium/long per mouse (WT vs APP) ----
    if ~isempty(COUNTS_BY_GENO)
        figure('Color','w'); hold on;

        types = {'Short','Medium','Long'};
        x = 1:numel(types);
        barWidth = 0.35;

        hasWT  = any(COUNTS_BY_GENO.genotype == "WT");
        hasAPP = any(COUNTS_BY_GENO.genotype == "APP");

        COL_WT  = [0.6 0.6 0.6];
        COL_APP = [0.39 0.58 0.93];

        meanWT  = nan(1,numel(types)); semWT  = nan(1,numel(types));
        meanAPP = nan(1,numel(types)); semAPP = nan(1,numel(types));

        for t = 1:numel(types)
            Tname  = types{t};
            colM   = ['mean' Tname];
            colS   = ['sem'  Tname];

            if hasWT
                rowWT = find(COUNTS_BY_GENO.genotype == "WT", 1);
                meanWT(t) = COUNTS_BY_GENO.(colM)(rowWT);
                semWT(t)  = COUNTS_BY_GENO.(colS)(rowWT);
            end
            if hasAPP
                rowAPP = find(COUNTS_BY_GENO.genotype == "APP", 1);
                meanAPP(t) = COUNTS_BY_GENO.(colM)(rowAPP);
                semAPP(t)  = COUNTS_BY_GENO.(colS)(rowAPP);
            end
        end

        if hasWT
            bar(x - barWidth/2, meanWT, barWidth, ...
                'FaceColor', COL_WT, 'EdgeColor','none');
            errorbar(x - barWidth/2, meanWT, semWT, 'k','LineStyle','none','LineWidth',1);
        end
        if hasAPP
            bar(x + barWidth/2, meanAPP, barWidth, ...
                'FaceColor', COL_APP, 'EdgeColor','none');
            errorbar(x + barWidth/2, meanAPP, semAPP, 'k','LineStyle','none','LineWidth',1);
        end

        xticks(x);
        xticklabels(types);
        xlabel('REM type');
        ylabel('# REM bouts per mouse (baseline)');
        title('Baseline REM bouts: Short / Medium / Long (WT vs APP)');
        if hasWT && hasAPP
            legend({'WT','APP'},'Location','northoutside','Orientation','horizontal');
        elseif hasWT
            legend({'WT'});
        elseif hasAPP
            legend({'APP'});
        end
        set(gca,'Box','off','FontSize',12);

        rem_counts_png = fullfile(S.out_dir, 'baseline_REM_counts_by_type_APPvsWT.png');
        saveas(gcf, rem_counts_png);
    end
end

% ===========================================
% 5) REM CLUSTER STRUCTURE (baseline only)
% ===========================================
CLUSTERS = table();
LONG_REM = table();
SUMMARY  = table();

if ~isempty(BASE_REM)
    rec_keys = unique(BASE_REM(:, {'file','mouse','genotype','condition'}));
    cluster_gap = S.cluster_gap_sec;

    for r = 1:height(rec_keys)
        fk = rec_keys.file(r);
        mk = rec_keys.mouse(r);
        gk = rec_keys.genotype(r);
        ck = rec_keys.condition(r);

        maskR = BASE_REM.file == fk & BASE_REM.mouse == mk;
        R = BASE_REM(maskR,:);
        if height(R) < 1, continue; end

        [~, idx] = sort(R.start_s);
        R = R(idx,:);

        nB = height(R);
        cluster_id = zeros(nB,1);
        cid = 1;
        cluster_id(1) = cid;
        for k = 2:nB
            gap = R.start_s(k) - R.end_s(k-1);
            if gap <= cluster_gap
                cluster_id(k) = cid;
            else
                cid = cid + 1;
                cluster_id(k) = cid;
            end
        end
        R.cluster_id = cluster_id;

        cIDs = unique(cluster_id);
        for cID = cIDs.'
            maskC = (R.cluster_id == cID);
            Rc = R(maskC,:);

            hasShort = any(Rc.rem_type == "Short");
            hasLong  = any(Rc.rem_type == "Long");

            if hasLong && hasShort
                ctype = "ShortLong";
            elseif hasLong && ~hasShort
                ctype = "LongOnly";
            elseif ~hasLong && hasShort
                ctype = "ShortOnly";
            else
                ctype = "MediumOnly"; % only Medium REMs
            end

            long_idx = find(Rc.rem_type == "Long");
            for li = reshape(long_idx,1,[])
                thisStart = Rc.start_s(li);
                n_short_before = nnz(Rc.rem_type == "Short" & Rc.start_s < thisStart);

                LONG_REM = [LONG_REM; ...
                    table(fk, mk, gk, ck, cID, ...
                          Rc.start_s(li), Rc.end_s(li), Rc.dur_s(li), ...
                          n_short_before, ctype, ...
                          'VariableNames', {'file','mouse','genotype','condition', ...
                                            'cluster_id','long_start_s','long_end_s','long_dur_s', ...
                                            'n_short_before_long','cluster_type'})]; %#ok<AGROW>
            end

            CLUSTERS = [CLUSTERS; ...
                table(fk, mk, gk, ck, cID, ctype, height(Rc), ...
                      'VariableNames', {'file','mouse','genotype','condition', ...
                                        'cluster_id','cluster_type','n_REM_in_cluster'})]; %#ok<AGROW>
        end
    end

    if ~isempty(CLUSTERS)
        mice2 = unique(CLUSTERS.mouse);
        for m = 1:numel(mice2)
            mID = mice2(m);
            maskM = CLUSTERS.mouse == mID;
            CM = CLUSTERS(maskM,:);

            geno = unique(CM.genotype);
            geno = geno(1);

            n_shortonly = nnz(CM.cluster_type == "ShortOnly");
            n_shortlong = nnz(CM.cluster_type == "ShortLong");
            n_longonly  = nnz(CM.cluster_type == "LongOnly");
            n_medonly   = nnz(CM.cluster_type == "MediumOnly");
            nC = height(CM);

            maskL = LONG_REM.mouse == mID;
            if any(maskL)
                mean_n_short_before_long = mean(LONG_REM.n_short_before_long(maskL),'omitnan');
            else
                mean_n_short_before_long = NaN;
            end

            SUMMARY = [SUMMARY; ...
                table(mID, geno, nC, ...
                      n_shortonly, n_shortlong, n_longonly, n_medonly, ...
                      mean_n_short_before_long, ...
                      'VariableNames', {'mouse','genotype','n_clusters', ...
                                        'n_shortonly','n_shortlong','n_longonly','n_mediumonly', ...
                                        'mean_n_short_before_long'})]; %#ok<AGROW>
        end
    end
end

clusters_csv = fullfile(S.out_dir, 'baseline_REM_clusters.csv');
longrem_csv  = fullfile(S.out_dir, 'baseline_long_REM_with_nShortBefore.csv');
summary_csv  = fullfile(S.out_dir, 'baseline_REM_cluster_summary_per_mouse.csv');
cluster_png  = '';

if ~isempty(CLUSTERS), writetable(CLUSTERS, clusters_csv); end
if ~isempty(LONG_REM),  writetable(LONG_REM,  longrem_csv); end
if ~isempty(SUMMARY),   writetable(SUMMARY,   summary_csv); end

if ~isempty(SUMMARY)
    propTbl = SUMMARY;
    propTbl.p_shortonly = propTbl.n_shortonly ./ propTbl.n_clusters;
    propTbl.p_shortlong = propTbl.n_shortlong ./ propTbl.n_clusters;
    propTbl.p_longonly  = propTbl.n_longonly  ./ propTbl.n_clusters;
    propTbl.p_medonly   = propTbl.n_mediumonly ./ propTbl.n_clusters;

    genos = unique(propTbl.genotype);
    types = {'shortonly','shortlong','longonly','medonly'};
    labels = {'ShortOnly','ShortLong','LongOnly','MediumOnly'};

    figure('Color','w'); hold on;
    x = 1:numel(types);
    barWidth = 0.35;
    COL_WT  = [0.6 0.6 0.6];
    COL_APP = [0.39 0.58 0.93];

    meanWT = nan(1,numel(types)); semWT = nan(1,numel(types));
    meanAPP = nan(1,numel(types)); semAPP = nan(1,numel(types));

    hasWT  = any(propTbl.genotype == "WT");
    hasAPP = any(propTbl.genotype == "APP");

    for t = 1:numel(types)
        pname = ['p_' types{t}];

        if hasWT
            maskWT = propTbl.genotype == "WT";
            vWT = propTbl.(pname)(maskWT);
            vWT = vWT(~isnan(vWT));
            if ~isempty(vWT)
                meanWT(t) = mean(vWT);
                semWT(t)  = std(vWT)/sqrt(numel(vWT));
            end
        end
        if hasAPP
            maskAPP = propTbl.genotype == "APP";
            vAPP = propTbl.(pname)(maskAPP);
            vAPP = vAPP(~isnan(vAPP));
            if ~isempty(vAPP)
                meanAPP(t) = mean(vAPP);
                semAPP(t)  = std(vAPP)/sqrt(numel(vAPP));
            end
        end
    end

    if hasWT
        bar(x - barWidth/2, meanWT, barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
        errorbar(x - barWidth/2, meanWT, semWT, 'k','LineStyle','none','LineWidth',1);
    end
    if hasAPP
        bar(x + barWidth/2, meanAPP, barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
        errorbar(x + barWidth/2, meanAPP, semAPP, 'k','LineStyle','none','LineWidth',1);
    end

    xticks(x);
    xticklabels(labels);
    ylabel('Proportion of REM clusters');
    xlabel('Cluster type');
    title('REM cluster structure (baseline): WT vs APP');
    if hasWT && hasAPP
        legend({'WT','APP'},'Location','northoutside','Orientation','horizontal');
    elseif hasWT
        legend({'WT'});
    elseif hasAPP
        legend({'APP'});
    end
    set(gca,'Box','off','FontSize',12);

    cluster_png = fullfile(S.out_dir, 'baseline_REM_cluster_types_APPvsWT.png');
    saveas(gcf, cluster_png);
end

% ======================
% 6) PACK OUTPUT
% ======================
OUT = struct();
OUT.ALL_REM          = ALL_REM;
OUT.BASE_REM         = BASE_REM;
OUT.REM_THRESH       = REM_THRESH;
OUT.COUNTS_PER_MOUSE = COUNTS_PER_MOUSE;
OUT.COUNTS_BY_GENO   = COUNTS_BY_GENO;
OUT.CLUSTERS         = CLUSTERS;
OUT.LONG_REM         = LONG_REM;
OUT.SUMMARY          = SUMMARY;
OUT.params           = S;
OUT.files = struct( ...
    'ALL_REM_csv',              all_rem_csv, ...
    'counts_per_mouse_csv',     counts_csv, ...
    'counts_by_geno_csv',       counts_geno_csv, ...
    'rem_counts_png',           rem_counts_png, ...
    'clusters_csv',             clusters_csv, ...
    'long_rem_csv',             longrem_csv, ...
    'summary_csv',              summary_csv, ...
    'cluster_png',              cluster_png);

fprintf('✅ REM length + structure analysis finished. Outputs in: %s\n', S.out_dir);
end

% -------------------------------------------------------------------------
% Local helper: robust filename parser for your naming scheme
% -------------------------------------------------------------------------
function info = parse_info_from_filename_relaxed(fname)
    % Expected shapes (examples):
    %   20251001-baseline-mouse1-APP_scores_1Hz.csv
    %   20251001_baseline_mouse1_APP_scored_scores_1Hz.csv
    %
    % We assume first 4 tokens (split by '-' or '_') are:
    %   date, condition, mouseID, genotype

    info = struct('date','','condition','','mouse','','genotype','','ok',false);

    [~, name, ~] = fileparts(fname);   % strip extension
    parts = regexp(name,'[-_]','split');

    if numel(parts) < 4
        return;
    end

    info.date      = parts{1};
    info.condition = lower(parts{2});
    info.mouse     = parts{3};
    info.genotype  = upper(parts{4});
    info.ok        = true;
end
