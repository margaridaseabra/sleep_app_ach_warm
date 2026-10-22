%% RUN_00_PREPARE_DATA
% Validate metadata and input files before analysis.

clear;
clc;

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);

addpath(fullfile(repoDir, 'config'));

P = paths_config();

fprintf('\n========================================\n');
fprintf('AMBIENT TEMPERATURE PROJECT\n');
fprintf('DATA CHECK\n');
fprintf('========================================\n\n');

if ~isfile(P.mice)
    error('Mouse metadata not found: %s', P.mice);
end

if ~isfile(P.sessions)
    error('Session metadata not found: %s', P.sessions);
end

Mice = readtable(P.mice, ...
    'TextType', 'string', ...
    'VariableNamingRule', 'preserve');

Sessions = readtable(P.sessions, ...
    'TextType', 'string', ...
    'VariableNamingRule', 'preserve');

fprintf('Mice: %d\n', height(Mice));
fprintf('Sessions: %d\n\n', height(Sessions));

missingMouseMetadata = setdiff( ...
    unique(Sessions.mouse_id), ...
    unique(Mice.mouse_id));

if ~isempty(missingMouseMetadata)
    fprintf('Sessions with missing mouse metadata:\n');
    disp(missingMouseMetadata);
    error('Some session mice are missing from config/mice.csv');
end

fprintf('Session counts:\n\n');

cohorts = ["thesis","new"];
conditions = ["baseline","warming","donepezil"];

for c = 1:numel(cohorts)
    for k = 1:numel(conditions)

        n = sum( ...
            Sessions.cohort == cohorts(c) & ...
            Sessions.condition == conditions(k));

        fprintf('%-8s %-10s : %d\n', ...
            cohorts(c), conditions(k), n);
    end
end

fprintf('\nChecking files...\n');

nMissing = 0;

for i = 1:height(Sessions)

    filePath = fullfile( ...
        P.data_root, ...
        char(Sessions.score_file(i)));

    if ~isfile(filePath)
        fprintf('MISSING: %s\n', filePath);
        nMissing = nMissing + 1;
    end
end

if nMissing > 0
    error('%d input scoring files are missing.', nMissing);
end

fprintf('All %d scoring files found.\n', height(Sessions));

fprintf('\nGenotype counts:\n');

[G, cohortGroup, genotypeGroup] = findgroups(Mice.cohort, Mice.genotype);
nMice = splitapply(@numel, Mice.mouse_id, G);

genotypeSummary = table( ...
    cohortGroup, ...
    genotypeGroup, ...
    nMice, ...
    'VariableNames', {'cohort','genotype','n_mice'});

disp(genotypeSummary);

fprintf('\nData check complete.\n');
