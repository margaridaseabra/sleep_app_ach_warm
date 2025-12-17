function ST = run_rm_anova_bouts(rows_perhr, state_str)
% run_rm_anova_bouts  Two-way RM ANOVA for bouts/hour in one state
%
%   ST = run_rm_anova_bouts(rows_perhr, 'REM')
%
% INPUT
%   rows_perhr : table from run_group_sleep_architecture (group_per_hour)
%   state_str  : 'WK','MA','NREM','REM'
%
% DESIGN
%   Within-subject factor: Time (hour_idx)
%   Between-subject factor: Genotype (WT vs APP)
%
% OUTPUT ST fields:
%   ST.state
%   ST.hours
%   ST.withinTbl        : ranova table
%   ST.betweenTbl       : anova(rm) table
%   ST.p_time
%   ST.p_interaction
%   ST.p_genotype

% ---------- 1) baseline & this state ----------
T = rows_perhr;
cond_str    = lower(strtrim(T.condition));
is_baseline = cond_str == "baseline";
T = T(is_baseline & T.state == state_str, :);

if isempty(T)
    error('No baseline rows for state %s.', state_str);
end

ST = struct();
ST.state = string(state_str);

% ---------- 2) Pivot to wide format ----------
all_hours = unique(T.hour_idx);
all_hours = sort(all_hours);
nH = numel(all_hours);

mice = unique(T.mouse);
nM   = numel(mice);

dataMat = nan(nM, nH);
genoVec = strings(nM,1);

for i = 1:nM
    mID = mice(i);
    maskMouse = (T.mouse == mID);

    g_this = unique(T.genotype(maskMouse));
    genoVec(i) = g_this(1);   % assume 1 genotype per mouse

    for h = 1:nH
        hr = all_hours(h);
        mask = maskMouse & (T.hour_idx == hr);
        vals = T.bouts_per_h(mask);
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            dataMat(i,h) = mean(vals);
        end
    end
end

hourNames = arrayfun(@(hr) sprintf('H%d', hr), all_hours, 'UniformOutput', false);
Twide = table(mice, genoVec, 'VariableNames', {'mouse','genotype'});
for h = 1:nH
    Twide.(hourNames{h}) = dataMat(:,h);
end

WithinDesign = table(all_hours(:), 'VariableNames', {'Time'});
measureStr   = sprintf('%s-%s', hourNames{1}, hourNames{end});

rm = fitrm(Twide, sprintf('%s ~ genotype', measureStr), ...
           'WithinDesign', WithinDesign);

withinTbl  = ranova(rm, 'WithinModel','Time');
betweenTbl = anova(rm);

ST.hours      = all_hours;
ST.withinTbl  = withinTbl;
ST.betweenTbl = betweenTbl;

% ---------- 3) Get term names robustly ----------
% withinTbl: term labels may be in a 'Term' column OR in RowNames
if any(strcmp(withinTbl.Properties.VariableNames,'Term'))
    wTerms = string(withinTbl.Term);
else
    wTerms = string(withinTbl.Properties.RowNames);
end

% p-column name: prefer pValueGG, else pValue, else last numeric column
if any(strcmp(withinTbl.Properties.VariableNames,'pValueGG'))
    pcol_within = 'pValueGG';
elseif any(strcmp(withinTbl.Properties.VariableNames,'pValue'))
    pcol_within = 'pValue';
else
    % crude fallback: last column
    pcol_within = withinTbl.Properties.VariableNames{end};
end

% betweenTbl: genotype term may be in Term or RowNames
if any(strcmp(betweenTbl.Properties.VariableNames,'Term'))
    bTerms = string(betweenTbl.Term);
else
    bTerms = string(betweenTbl.Properties.RowNames);
end

if any(strcmp(betweenTbl.Properties.VariableNames,'pValue'))
    pcol_between = 'pValue';
else
    pcol_between = betweenTbl.Properties.VariableNames{end};
end

% ---------- 4) Extract p-values ----------
idxTime = find(wTerms == "Time", 1);
idxInt  = find(wTerms == "Time:genotype", 1);

if ~isempty(idxTime)
    ST.p_time = withinTbl.(pcol_within)(idxTime);
else
    ST.p_time = NaN;
end

if ~isempty(idxInt)
    ST.p_interaction = withinTbl.(pcol_within)(idxInt);
else
    ST.p_interaction = NaN;
end

idxGen = find(bTerms == "genotype", 1);
if ~isempty(idxGen)
    ST.p_genotype = betweenTbl.(pcol_between)(idxGen);
else
    ST.p_genotype = NaN;
end

% ---------- 5) Print summary ----------
fprintf('\n===== Two-way RM ANOVA (Time × Genotype) for %s bouts/hour =====\n', state_str);
if ~isnan(ST.p_time)
    fprintf('  Time effect:           p = %.4g\n', ST.p_time);
else
    fprintf('  Time effect:           p = NaN (could not extract)\n');
end
if ~isnan(ST.p_interaction)
    fprintf('  Time × Genotype:       p = %.4g\n', ST.p_interaction);
else
    fprintf('  Time × Genotype:       p = NaN (could not extract)\n');
end
if ~isnan(ST.p_genotype)
    fprintf('  Genotype main effect:  p = %.4g\n', ST.p_genotype);
else
    fprintf('  Genotype main effect:  p = NaN (could not extract)\n');
end
end
