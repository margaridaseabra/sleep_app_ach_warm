function OUT = run_group_sleep_architecture(input_dir, varargin)
% RUN_GROUP_SLEEP_ARCHITECTURE
% Batch sleep-architecture analysis + genotype/condition comparisons.
%
% Assumes filenames like:
%   <date>-<condition>-<mouse>-<APP_or_WT>_scores_1Hz.csv
% Example:
%   20251001-baseline-mouse1-APP_scores_1Hz.csv
%
% Usage:
% OUT = run_group_sleep_architecture('/path/to/scores', ...
%   'pattern','*_scores_1Hz.csv', ...
%   'codes', struct('WK',0,'NREM',1,'REM',2,'MA',15), ...
%   'includeMA', false, ...
%   'ma_thresh_sec', 15, ...
%   'reclassify_short_wake_to_MA', true, ...
%   'out_dir', fullfile(pwd,'group_out'));

% ----------------- args -----------------
p = inputParser;
addRequired(p,'input_dir',@ischar);
addParameter(p,'pattern','*_scores_1Hz.csv',@ischar);
addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
addParameter(p,'includeMA',false,@islogical);
addParameter(p,'ma_thresh_sec',15,@(x)isscalar(x)&&x>0);
addParameter(p,'reclassify_short_wake_to_MA',true,@islogical);
addParameter(p,'out_dir','',@ischar);
parse(p, input_dir, varargin{:});
S = p.Results;

% ---- guard: function must be on path
if exist('sleep_architecture_from_scores','file') ~= 2
    error(['sleep_architecture_from_scores.m is not on the MATLAB path.\n' ...
           'Add it, e.g.: addpath(genpath(''%s''))'], pwd);
end

assert(isfolder(S.input_dir), 'Input folder not found: %s', S.input_dir);
if isempty(S.out_dir), S.out_dir = fullfile(S.input_dir, 'group_out'); end
if ~isfolder(S.out_dir), mkdir(S.out_dir); end

% ----------------- discover files -----------------
F = dir(fullfile(S.input_dir, '**', S.pattern));
F = F(~[F.isdir]);
assert(~isempty(F), 'No files matched "%s" under %s', S.pattern, S.input_dir);

% ----------------- fixed collector schemas -----------------
schema_overall_names = {'state','n_bouts','total_dur_s','mean_bout_dur_s', ...
                        'file','date','condition','mouse','genotype'};
schema_overall_types = {'string','double','double','double','string','string','string','string','string'};
rows_overall = table('Size',[0 numel(schema_overall_names)], ...
    'VariableTypes', schema_overall_types, ...
    'VariableNames', schema_overall_names);

schema_perhr_names = {'hour_idx','hour_start_s','dur_s','bouts_per_h','mean_bout_dur_s', ...
                      'state','file','date','condition','mouse','genotype'};
schema_perhr_types = {'double','double','double','double','double', ...
                      'string','string','string','string','string','string'};
rows_perhr = table('Size',[0 numel(schema_perhr_names)], ...
    'VariableTypes', schema_perhr_types, ...
    'VariableNames', schema_perhr_names);

schema_meta_names = {'file','date','condition','mouse','genotype','ok'};
schema_meta_types = {'string','string','string','string','string','logical'};
META = table('Size',[0 numel(schema_meta_names)], ...
    'VariableTypes', schema_meta_types, ...
    'VariableNames', schema_meta_names);

% ----------------- loop files -----------------
for i = 1:numel(F)
    csv_path = fullfile(F(i).folder, F(i).name);
    info = parse_info_from_filename(F(i).name);

    % always track metadata (standardized columns only)
    META = [META; {string(F(i).name), string(info.date), string(info.condition), ...
                   string(info.mouse), string(info.genotype), logical(info.ok)}]; %#ok<AGROW>

    if ~info.ok
        warning('Filename not recognized pattern, skipping: %s', F(i).name);
        continue;
    end

    base_prefix = erase(F(i).name, '_scores_1Hz.csv'); % for per-file exports

    % ---- run per-file analysis
    try
        OUT_i = sleep_architecture_from_scores(csv_path, ...
            'codes', S.codes, ...
            'includeMA', S.includeMA, ...
            'ma_thresh_sec', S.ma_thresh_sec, ...
            'reclassify_short_wake_to_MA', S.reclassify_short_wake_to_MA, ...
            'out_prefix', base_prefix, ...
            'out_dir', S.out_dir, ...
            'write_log', true, ...
            'verbose', false);
    catch ME
        warning('Failed file %s: %s', F(i).name, ME.message);
        continue;
    end

    % ---- OVERALL (normalized row build)
    O = OUT_i.overall;
    if istable(O) && ~isempty(O)
        O2 = table( ...
            string(O.state), ...
            double(O.n_bouts), ...
            double(O.total_dur_s), ...
            double(O.mean_bout_dur_s), ...
            repmat(string(F(i).name),     height(O),1), ...
            repmat(string(info.date),     height(O),1), ...
            repmat(string(info.condition),height(O),1), ...
            repmat(string(info.mouse),    height(O),1), ...
            repmat(string(info.genotype), height(O),1), ...
            'VariableNames', schema_overall_names);
        rows_overall = [rows_overall; O2]; %#ok<AGROW>
    end

    % ---- PER-HOUR (normalized row build)
    P = OUT_i.per_hour;
    if istable(P) && ~isempty(P)
        states = {'wk','nrem','rem'};
        if S.includeMA && any(strcmpi(OUT_i.overall.state,'MA')), states = [states, {'ma'}]; end

        for s = 1:numel(states)
            st = states{s};
            dur_col = [st '_dur_s'];
            bph_col = [st '_bouts_per_h'];
            md_col  = [st '_mean_bout_dur_s'];
            if ~all(ismember({dur_col,bph_col}, P.Properties.VariableNames)), continue; end

            mdur = NaN(height(P),1);
            if ismember(md_col, P.Properties.VariableNames)
                mdur = double(P.(md_col));
            end

            Trow2 = table( ...
                double(P.hour_idx), ...
                double(P.hour_start_s), ...
                double(P.(dur_col)), ...
                double(P.(bph_col)), ...
                mdur, ...
                repmat(string(upper(st)),        height(P),1), ...
                repmat(string(F(i).name),        height(P),1), ...
                repmat(string(info.date),        height(P),1), ...
                repmat(string(info.condition),   height(P),1), ...
                repmat(string(info.mouse),       height(P),1), ...
                repmat(string(info.genotype),    height(P),1), ...
                'VariableNames', schema_perhr_names);

            rows_perhr = [rows_perhr; Trow2]; %#ok<AGROW>
        end
    end
end

% remove any per-hour summary rows if they exist (hour_idx == -1)
if ~isempty(rows_perhr) && any(rows_perhr.hour_idx == -1)
    rows_perhr(rows_perhr.hour_idx == -1,:) = [];
end

% ----------------- save group tables -----------------
overall_csv = fullfile(S.out_dir, 'group_overall.csv');
perhour_csv = fullfile(S.out_dir, 'group_per_hour.csv');
meta_csv    = fullfile(S.out_dir, 'group_meta.csv');

if ~isempty(rows_overall), writetable(rows_overall, overall_csv); else
    warning('Overall table empty; not writing CSV.');
end
if ~isempty(rows_perhr), writetable(rows_perhr, perhour_csv); else
    warning('Per-hour table empty; not writing CSV.');
end
if ~isempty(META), writetable(META, meta_csv); end

% ----------------- plots (only if we have data) -----------------
if ~isempty(rows_overall) && ~isempty(rows_perhr)
    make_group_plots(rows_overall, rows_perhr, S.out_dir);
end

% ----------------- return -----------------
OUT = struct('overall',rows_overall,'per_hour',rows_perhr,'meta',META, ...
             'out_dir',S.out_dir,'files',struct('overall_csv',overall_csv, ...
                                                'perhour_csv',perhour_csv, ...
                                                'meta_csv',meta_csv));
fprintf('✅ Group run finished. Outputs in: %s\n', S.out_dir);
end

% ===================== helpers =====================

function info = parse_info_from_filename(fname)
% Expect: <date>-<condition>-<mouse>-<APP_or_WT>_scores_1Hz.csv
% Returns only standardized fields used by META (no "raw" field).
info = struct('ok',false,'date','','condition','','mouse','','genotype','');
name = erase(fname, '_scores_1Hz.csv');
parts = split(name, '-');
if numel(parts) < 4, return; end
info.date      = string(parts{1});
info.condition = string(parts{2});
info.mouse     = string(parts{3});
geno = upper(string(parts{4}));
if contains(geno,'APP'), info.genotype = "APP";
elseif contains(geno,'WT'), info.genotype = "WT";
else, info.genotype = geno;
end
info.ok = true;
end

function make_group_plots(OVERALL, PERHOUR, out_dir)
% Simple, consistent palette per state
COL = struct('WK',[0.55 0.55 0.55], 'NREM',[0.30 0.50 0.85], 'REM',[0.80 0.35 0.35], 'MA',[0.60 0.60 0.20]);

% which states exist?
states = unique(OVERALL.state,'stable'); states = cellstr(states)';

% -------- Fig 1: Total duration (min) by GENOTYPE x CONDITION
for s = 1:numel(states)
    st = states{s};
    sub = OVERALL(strcmp(OVERALL.state, st), :);
    if isempty(sub), continue; end
    sub.total_min = sub.total_dur_s/60;
    G = agg_mean_sem(sub, {'genotype','condition'}, 'total_min');
    if isempty(G), continue; end
    f = figure('Color','w','Name', ['Total ' st ' (min)']); hold on
    cats = strcat(G.genotype, " | ", G.condition);
    b = bar(categorical(cats), G.mean); set(b,'FaceColor', pick_col(COL, st));
    errorbar(categorical(cats), G.mean, G.sem, 'k.', 'LineWidth',1);
    ylabel('Total duration (min)'); title(sprintf('Total %s (min)', st));
    xtickangle(30); box off
    saveas(f, fullfile(out_dir, sprintf('fig_total_%s.png', lower(st))));
end

% -------- Fig 2: Number of bouts by GENOTYPE x CONDITION
for s = 1:numel(states)
    st = states{s};
    sub = OVERALL(strcmp(OVERALL.state, st), :);
    if isempty(sub), continue; end
    G = agg_mean_sem(sub, {'genotype','condition'}, 'n_bouts');
    if isempty(G), continue; end
    f = figure('Color','w','Name', ['Bouts ' st]); hold on
    cats = strcat(G.genotype, " | ", G.condition);
    b = bar(categorical(cats), G.mean); set(b,'FaceColor', pick_col(COL, st));
    errorbar(categorical(cats), G.mean, G.sem, 'k.', 'LineWidth',1);
    ylabel('Number of bouts'); title(sprintf('Bouts — %s', st));
    xtickangle(30); box off
    saveas(f, fullfile(out_dir, sprintf('fig_bouts_%s.png', lower(st))));
end

% -------- Fig 3: Mean bout duration (s) by GENOTYPE x CONDITION
for s = 1:numel(states)
    st = states{s};
    sub = OVERALL(strcmp(OVERALL.state, st), :);
    if isempty(sub), continue; end
    G = agg_mean_sem(sub, {'genotype','condition'}, 'mean_bout_dur_s');
    if isempty(G), continue; end
    f = figure('Color','w','Name', ['Mean bout dur ' st]); hold on
    cats = strcat(G.genotype, " | ", G.condition);
    b = bar(categorical(cats), G.mean); set(b,'FaceColor', pick_col(COL, st));
    errorbar(categorical(cats), G.mean, G.sem, 'k.', 'LineWidth',1);
    ylabel('Mean bout duration (s)'); title(sprintf('Mean bout duration — %s', st));
    xtickangle(30); box off
    saveas(f, fullfile(out_dir, sprintf('fig_mean_bout_%s.png', lower(st))));
end

% -------- Fig 4: Bouts per hour trajectories (mean ± SEM across files)
PH = PERHOUR; PH = PH(~isnan(PH.bouts_per_h),:);
if ~isempty(PH)
    % aggregate mean/sem by (state, genotype, condition, hour_idx)
    Tagg = agg_mean_sem(PH, {'state','genotype','condition','hour_idx'}, 'bouts_per_h');
    Tagg.Properties.VariableNames{'hour_idx'} = 'hour';  % nicer name

    for s = 1:numel(states)
        st = states{s};
        sub = Tagg(strcmp(Tagg.state, st), :);
        if isempty(sub), continue; end

        f = figure('Color','w','Name', ['Bouts per hour — ' st]); 
        tiledlayout('flow');

        Ugeno = unique(sub.genotype,'stable'); 
        Ucond = unique(sub.condition,'stable');

        for g = 1:numel(Ugeno)
            for c = 1:numel(Ucond)
                sc = sub(strcmp(sub.genotype,Ugeno{g}) & strcmp(sub.condition,Ucond{c}), :);
                if isempty(sc), continue; end
                nexttile; hold on
                errorbar(double(sc.hour), sc.mean, sc.sem, 'o-','LineWidth',1.25, ...
                         'Color', pick_col(COL, st));
                title(sprintf('%s | %s', Ugeno{g}, Ucond{c}));
                xlabel('Hour'); ylabel('Bouts per hour');
                grid on; box off
            end
        end

        sgtitle(sprintf('Bouts per hour — %s', st));
        saveas(f, fullfile(out_dir, sprintf('fig_bph_byhour_%s.png', lower(st))));
    end
end

end % Missing end added here

% function G = agg_mean_sem(T, groupVars, valueVar)
% % Robust mean ± SEM aggregator using findgroups + splitapply
% % Returns a table G that includes the grouping columns + mean + sem.
% 
% % Pull the value vector
% vals = double(T.(valueVar));
% 
% % Build grouping vectors
% gv = cell(1, numel(groupVars));
% for k = 1:numel(groupVars)
%     v = T.(groupVars{k});
%     % normalize text columns to string for stable grouping
%     if ischar(v) || isstring(v) || iscellstr(v)
%         gv{k} = string(v);
%     else
%         gv{k} = v;
%     end
% end
% 
% 
% % Group + aggregate
% [Gid, keys{:}] = findgroups(gv{:}); %#ok<AGROW,CCAT1>
% m   = splitapply(@(x) mean(x,'omitnan'), vals, Gid);
% sd  = splitapply(@(x) std(x,'omitnan'),  vals, Gid);
% n   = splitapply(@(x) sum(~isnan(x)),    vals, Gid);
% sem = sd ./ max(sqrt(n), 1);
% 
% % Build output table with key columns preserved
% G = table();
% for k = 1:numel(groupVars)
%     G.(groupVars{k}) = keys{k};
% end
% G.mean = m;
% G.sem  = sem;
% end

function c = pick_col(COL, st)
switch upper(st)
    case 'WK',   c = COL.WK;
    case 'NREM', c = COL.NREM;
    case 'REM',  c = COL.REM;
    case 'MA',   c = COL.MA;
    otherwise,   c = [0.5 0.5 0.5];
end
end

function G = agg_mean_sem(T, groupVars, valueVar)
% Robust mean ± SEM aggregator using findgroups + splitapply
% Returns a table G with columns: groupVars..., mean, sem

% Values
vals = double(T.(valueVar));

% Normalize grouping vectors
gv = cell(1, numel(groupVars));
for k = 1:numel(groupVars)
    v = T.(groupVars{k});
    if ischar(v) || isstring(v) || iscellstr(v)
        gv{k} = string(v);
    else
        gv{k} = v;
    end
end

% Preallocate receivers for findgroups varargout
keysOut = cell(1, numel(groupVars));
[Gid, keysOut{:}] = findgroups(gv{:});

% Aggregate
m  = splitapply(@(x) mean(x,'omitnan'), vals, Gid);
sd = splitapply(@(x) std(x,'omitnan'),  vals, Gid);
n  = splitapply(@(x) sum(~isnan(x)),    vals, Gid);
sem = sd ./ max(sqrt(n), 1);

% Build output table
G = table();
for k = 1:numel(groupVars)
    G.(groupVars{k}) = keysOut{k};
end
G.mean = m;
G.sem  = sem;
end
