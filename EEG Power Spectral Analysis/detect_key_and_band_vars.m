function [mouseVar, genoVar, condVar, stateVar, bandVars] = detect_key_and_band_vars(T)

vars = T.Properties.VariableNames;

% Mouse
candidates_mouse = {'mouse','mouseID','Mouse','mouse_id'};
mouseVar = find_first_existing(vars, candidates_mouse);

% Genotype
candidates_geno = {'geno','genotype','Genotype'};
genoVar  = find_first_existing(vars, candidates_geno);

% Condition
candidates_cond = {'cond','condition','Condition'};
condVar  = find_first_existing(vars, candidates_cond);

% State
candidates_state = {'state','State','stage','Stage'};
stateVar = find_first_existing(vars, candidates_state);

if any(cellfun(@isempty, {mouseVar,genoVar,condVar,stateVar}))
    error('Could not auto-detect key variables: mouse/genotype/condition/state.');
end

% Band variables = numeric variables that are not keys and not obvious index vars
isNum = varfun(@isnumeric, T, 'OutputFormat','uniform');
numVars = vars(isNum);

exclude = {mouseVar, genoVar, condVar, stateVar, ...
           'epoch','epoch_idx','t','time_s','sample_idx','fs'};
bandVars = setdiff(numVars, exclude);

if isempty(bandVars)
    error('No numeric band-power variables detected.');
end
end

function v = find_first_existing(allVars, candidates)
v = '';
for k = 1:numel(candidates)
    if any(strcmp(allVars, candidates{k}))
        v = candidates{k};
        return;
    end
end
end
