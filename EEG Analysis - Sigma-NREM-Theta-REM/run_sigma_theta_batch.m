function run_sigma_theta_batch(eegDir, scoreDir, outRoot)
% run_sigma_theta_batch(eegDir, scoreDir, outRoot)
%
% Runs psd_sigma_theta_analysis for all EEG/score pairs,
% groups results across mice, genotypes and conditions, and
% saves:
%   - per-session outputs under outRoot/genotype/condition/mouseID
%   - a group CSV: outRoot/SigmaTheta_modulation_allmice.csv
%
% Assumes:
%   EEG .mat names like:  20251023-drugs_mouse3_notched50Hz_bw3.mat
%   CSV scores like:      20251023-drugs-mouse3-WT_scored_scores_1Hz.csv
%
% Uses parse_meta_from_scorefile.m to read condition/mouse/genotype.

    if nargin < 3 || isempty(outRoot)
        outRoot = 'SigmaTheta_ModAnalysis';
    end
    if ~exist(outRoot,'dir'); mkdir(outRoot); end

    % ---- SLEEP CODES: EDIT if your numeric codes differ ----
    % Example: 1=Wake, 4=NREM, 9=REM, 15=MA
    codes = struct('WK',0,'NREM',1,'REM',2,'MA',15);

    % ---- Find EEG files ----
    eegFiles = dir(fullfile(eegDir, '*notched50Hz_bw3*.mat'));
    if isempty(eegFiles)
        error('No EEG files found in %s', eegDir);
    end

    allTbl = table();

    for k = 1:numel(eegFiles)
        eegFile = fullfile(eegDir, eegFiles(k).name);
        [~, eegName] = fileparts(eegFiles(k).name);

        % 1) Remove "_notched..." suffix
        eegBase = regexprep(eegName, '_notched.*$', '');
        %    e.g. 20251023-drugs_mouse3

        % 2) Replace "_" with "-" to match CSV style
        patternBase = strrep(eegBase, '_', '-');
        %    e.g. 20251023-drugs-mouse3

        % 3) Look for matching score CSV
        scorePattern = sprintf('%s-*_scored_scores_1Hz*.csv', patternBase);
        scoreCandidates = dir(fullfile(scoreDir, scorePattern));

        if isempty(scoreCandidates)
            warning('No score CSV found for EEG %s with pattern %s', ...
                    eegName, scorePattern);
            continue;
        elseif numel(scoreCandidates) > 1
            warning('Multiple score CSVs for EEG %s, using first.', eegName);
        end

        scoreFile = fullfile(scoreDir, scoreCandidates(1).name);

        % ---- Parse mouseID, condition, genotype from CSV name ----
        meta = parse_meta_from_scorefile(scoreFile);
        mouseID   = meta.mouseID;
        condition = meta.condition;   % e.g. baseline / ambtemp / drugs
        genotype  = meta.genotype;    % WT / APP

        % ---- Session-specific output directory ----
        sesOutDir = fullfile(outRoot, genotype, condition, mouseID);
        if ~exist(sesOutDir,'dir'); mkdir(sesOutDir); end

        outPrefix = sprintf('%s_%s_%s', mouseID, condition, genotype);

        fprintf('\n=== Sigma/Theta modulation: %s | %s | %s ===\n', ...
                mouseID, genotype, condition);

        % ---- Run the main analysis for this session ----
        OUT = psd_sigma_theta_analysis( ...
            eegFile, scoreFile, ...
            'codes', codes, ...
            'mouse_id',  mouseID, ...
            'session',   char(condition), ...
            'out_prefix', outPrefix, ...
            'out_dir',   sesOutDir, ...
            'verbose',   true);

        % Save full OUT struct
        save(fullfile(sesOutDir, 'SigmaTheta_OUT.mat'), 'OUT');

        % ---- Append one row to the group summary table ----
        newRow = table( ...
            string(mouseID), ...
            string(genotype), ...
            string(condition), ...
            OUT.sigma.mod_power, ...
            OUT.sigma.mod_peak_f, ...
            OUT.sigma.mod_peak_amp, ...
            OUT.theta.mod_power, ...
            OUT.theta.mod_peak_f, ...
            OUT.theta.mod_peak_amp, ...
            'VariableNames', {'MouseID','Genotype','Condition', ...
                              'sigma_mod_power','sigma_mod_peak_f','sigma_mod_peak_amp', ...
                              'theta_mod_power','theta_mod_peak_f','theta_mod_peak_amp'});

        allTbl = [allTbl; newRow]; %#ok<AGROW>
    end

    % ---- Save big CSV for stats ----
    outCsv = fullfile(outRoot, 'SigmaTheta_modulation_allmice.csv');
    writetable(allTbl, outCsv);
    fprintf('\n📁 Saved group sigma/theta modulation summary to:\n   %s\n', outCsv);
end
