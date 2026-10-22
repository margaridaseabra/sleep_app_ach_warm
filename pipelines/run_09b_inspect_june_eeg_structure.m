%% run_09b_inspect_june_eeg_structure.m
% Inspect the internal structure of one June MAT file without changing data.
%
% Run after run_09_validate_eeg_sources.m.
%
% It prints:
%   - top-level variable names
%   - class and dimensions
%   - fields of top-level structs/cells when safe to inspect
%
% No files are modified.

clear;
clc;

dataRoot = fullfile(getenv('HOME'),'Desktop','AmbTemp_proj');
f = fullfile(dataRoot,'scored_files_june26batch', ...
    '290526-m59-bas-1_final_scoring_dash.mat');

assert(isfile(f),'File not found:\n%s',f);

fprintf('\n============================================================\n');
fprintf('JUNE MAT STRUCTURE INSPECTION\n');
fprintf('============================================================\n');
fprintf('%s\n\n',f);

info = whos('-file',f);

fprintf('TOP-LEVEL VARIABLES\n');
fprintf('------------------------------------------------------------\n');
for i = 1:numel(info)
    fprintf('%-35s  %-12s  %s  %.2f MB\n', ...
        info(i).name, ...
        info(i).class, ...
        mat2str(info(i).size), ...
        info(i).bytes/1024/1024);
end

fprintf('\nSTRUCT / CELL CONTENTS\n');
fprintf('------------------------------------------------------------\n');

for i = 1:numel(info)
    nm = info(i).name;

    if strcmp(info(i).class,'struct')
        fprintf('\nSTRUCT: %s\n',nm);
        try
            tmp = load(f,nm);
            x = tmp.(nm);
            print_struct_fields(x,nm,0,3);
        catch ME
            fprintf('  Could not inspect: %s\n',ME.message);
        end

    elseif strcmp(info(i).class,'cell') && info(i).bytes < 50*1024*1024
        fprintf('\nCELL: %s\n',nm);
        try
            tmp = load(f,nm);
            x = tmp.(nm);
            fprintf('  size = %s\n',mat2str(size(x)));
            nshow = min(numel(x),10);
            for j = 1:nshow
                fprintf('  {%d}: class=%s size=%s\n', ...
                    j,class(x{j}),mat2str(size(x{j})));
            end
        catch ME
            fprintf('  Could not inspect: %s\n',ME.message);
        end
    end
end

fprintf('\nLIKELY EEG / SAMPLING-RATE NAMES\n');
fprintf('------------------------------------------------------------\n');

names = string({info.name});
tokens = ["eeg","EEG","signal","data","record","stream","fs","freq", ...
          "sampling","sample_rate","samplerate","frequency"];

for i = 1:numel(names)
    low = lower(names(i));
    if any(contains(low,lower(tokens)))
        fprintf('  %s\n',names(i));
    end
end

fprintf('\n============================================================\n');
fprintf('DONE\n');
fprintf('============================================================\n');

function print_struct_fields(x,prefix,depth,maxDepth)

if depth >= maxDepth || isempty(x) || ~isstruct(x)
    return;
end

f = fieldnames(x);

for k = 1:numel(f)
    name = f{k};

    try
        y = x(1).(name);
        fprintf('%s%s.%s : class=%s size=%s\n', ...
            repmat('  ',1,depth+1),prefix,name,class(y),mat2str(size(y)));

        if isstruct(y)
            print_struct_fields(y,[prefix '.' name],depth+1,maxDepth);
        end
    catch
        fprintf('%s%s.%s : <unable to inspect>\n', ...
            repmat('  ',1,depth+1),prefix,name);
    end
end
end
