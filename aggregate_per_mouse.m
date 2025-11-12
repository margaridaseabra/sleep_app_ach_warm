function M = aggregate_per_mouse(ALL)
% Per-mouse/condition aggregation (median of numeric features, omit NaNs).
% Avoids pseudo-replication by summarizing bouts within each mouse×condition.

% Grouping (works with string/categorical/cellstr)
[G, mice, conds] = findgroups(ALL.mouse, ALL.condition);

% Start output with the unique group keys
M = table(mice, conds, 'VariableNames', {'mouse','condition'});

% Pick numeric feature columns (exclude metadata if present)
meta = {'source_file','mouse','condition','bout'};
isNum = varfun(@isnumeric, ALL, 'OutputFormat','uniform');
numNames = setdiff(ALL.Properties.VariableNames(isNum), meta, 'stable');

% Aggregate each numeric feature with a groupwise median
for i = 1:numel(numNames)
    v = numNames{i};
    M.(v) = splitapply(@(x) median(x,'omitnan'), ALL.(v), G);
end
end
