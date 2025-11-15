function run_all_scored_files()
% Batch: scored .mat -> REM selection -> REM TFR analysis
%
% Assumptions match your current functions:
% - select_REM_episodes saves to:
%     N:\SUN-IN-Kjaerby\0_Personal folders\Margarida\REM analysis
%   with name: <base>_REMselection.mat   (or sometimes *_REMselection*.mat)
% - rem_timefreq_analysis reads that file and writes <base>_REM_TFR\*
%
% Edit only the two folder paths below if needed.

scored_dir = 'N:\SUN-IN-Kjaerby\0_Personal folders\Margarida\scored files';
sel_out_dir = 'N:\SUN-IN-Kjaerby\0_Personal folders\Margarida\REM analysis';  % matches your select_ function

method = 'spectrogram';
normopt = 'none';

files = dir(fullfile(scored_dir, '*.mat'));
fprintf('Found %d .mat files in "%s"\n', numel(files), scored_dir);

for i = 1:numel(files)
    in_file = fullfile(files(i).folder, files(i).name);

    % skip already-processed selection / analysis files inside the scored folder (just in case)
    if contains(files(i).name, '_REMselection', 'IgnoreCase', true) || ...
       contains(files(i).name, '_REM_TFR', 'IgnoreCase', true)
        fprintf('[%2d/%2d] Skip (already looks like output): %s\n', i, numel(files), files(i).name);
        continue;
    end

    fprintf('\n[%2d/%2d] Processing: %s\n', i, numel(files), files(i).name);

    try
        % ---- Step 1: REM selection ----
        fprintf(' -> select_REM_episodes...\n');
        select_REM_episodes('file', in_file);  % your function does not return a path

        % The REM selection file is written to sel_out_dir with base+"_REMselection*.mat"
        [~, base] = fileparts(in_file);
        sel_file = find_latest_selection(sel_out_dir, base);
        if sel_file == ""
            warning('   Could not locate a REMselection file for %s. Skipping analysis.', base);
            continue;
        else
            fprintf('    Found selection file: %s\n', sel_file);
        end

        % ---- Step 2: Time-frequency analysis ----
        fprintf(' -> rem_timefreq_analysis (%s, %s)...\n', method, normopt);
        rem_timefreq_analysis('file', sel_file, 'method', method, 'norm', normopt);

        fprintf('    [OK] Done: %s\n', base);
    catch ME
        warning('    [FAILED] %s\n    -> %s', files(i).name, ME.message);
    end
end

fprintf('\nAll done.\n');
end

function sel_file = find_latest_selection(sel_out_dir, base)
% Look for any "<base>_REMselection*.mat" under sel_out_dir and return the newest
d = dir(fullfile(sel_out_dir, base + "_REMselection*.mat"));
if isempty(d)
    sel_file = "";
else
    [~, idx] = max([d.datenum]);
    sel_file = string(fullfile(d(idx).folder, d(idx).name));
end
end
