function plot_group_psd_curves_outdated(baseOutDir)
% plot_group_psd_curves(baseOutDir)
%
% baseOutDir should be the same root you used before:
%   e.g. baseOutDir = 'EEG_PSD_AllMice';
%
% Assumes folder structure:
%   baseOutDir/genotype/condition/mouseID/PSD_data.mat
%
% For each condition, creates a figure with 3 subplots:
%   Wake, NREM, REM
% Each subplot shows WT vs APP mean PSD (in dB) with SEM shading.

    if nargin < 1 || isempty(baseOutDir)
        baseOutDir = 'EEG_PSD_AllMice';
    end

    % Genotype codes as they appear in folder names:
    genotypes = {'WT','APP'};   % change 'APP' to 'APPPS1' if that's your folder name

    % How we want to label them in the legend:
    labelMap = struct('WT','WT', 'APP','APP/PS1');

    % States and their names in the PSD_data struct:
    states = {'Wake','NREM','REM'};

    % Colors (roughly grey for WT, orange for APP)
    colMap.WT  = [0.3 0.3 0.3];
    colMap.APP = [1.0 0.6 0.1];

    % ---------- Collect all condition names ----------
    condSet = {};
    for g = 1:numel(genotypes)
        gDir = fullfile(baseOutDir, genotypes{g});
        if ~exist(gDir, 'dir'), continue; end
        d = dir(gDir);
        d = d([d.isdir] & ~ismember({d.name},{'.','..'}));
        condSet = [condSet, {d.name}]; %#ok<AGROW>
    end
    conditions = unique(condSet);

    if isempty(conditions)
        error('No condition folders found in %s', baseOutDir);
    end

    % ---------- Loop over conditions ----------
    for ci = 1:numel(conditions)
        cond = conditions{ci};

        figure('Color','w','Name',cond);
        nStates = numel(states);

        for si = 1:nStates
            st = states{si};
            subplot(1, nStates, si); hold on;

            hLines = gobjects(1, numel(genotypes));  % to store line handles

for gi = 1:numel(genotypes)
    geno = genotypes{gi};
    gDir  = fullfile(baseOutDir, geno, cond);
    if ~exist(gDir,'dir')
        warning('Missing folder: %s', gDir);
        continue;
    end

    mList = dir(gDir);
    mList = mList([mList.isdir] & ~ismember({mList.name},{'.','..'}));

    allPSD_dB = [];
    F_ref     = [];

    for mi = 1:numel(mList)
        sesDir  = fullfile(gDir, mList(mi).name);
        matFile = fullfile(sesDir, 'PSD_data.mat');
        if ~exist(matFile,'file'), continue; end

        S = load(matFile);
        if ~isfield(S,'result') || ~isfield(S.result, st), continue; end

        F   = S.result.(st).F;
        PSD = S.result.(st).PSD;

        PSDdB = 10*log10(PSD);

        if isempty(F_ref)
            F_ref = F;
        elseif numel(F) ~= numel(F_ref) || any(F ~= F_ref)
            warning('Frequency mismatch in %s, skipping.', matFile);
            continue;
        end

        allPSD_dB(:, end+1) = PSDdB; %#ok<AGROW>
    end

    if isempty(allPSD_dB)
        warning('No PSD data for %s, %s, %s', cond, st, geno);
        continue;
    end

    mPSD  = mean(allPSD_dB, 2, 'omitnan');
    nSubs = sum(~isnan(allPSD_dB(1,:)));
    sPSD  = std(allPSD_dB, 0, 2, 'omitnan') ./ max(nSubs,1).^0.5;

    c = colMap.(geno);

    x = F_ref;
    upper = mPSD + sPSD;
    lower = mPSD - sPSD;

    % Shaded SEM (very light, no edge)
    fill([x; flipud(x)], [upper; flipud(lower)], c, ...
         'FaceAlpha', 0.20, 'EdgeColor', 'none');

    % Mean line (thicker, slightly darker)
    hLines(gi) = plot(x, mPSD, 'Color', c*0.9, 'LineWidth', 2.0);
end

xlim([0 100]);
xlabel('Frequency (Hz)');
if si == 1
    ylabel('Power (dB)');
end
title(st);
grid on;

% Only once per condition, e.g. on the REM panel
if si == numel(states)
    % use the line handles, not the patches
    valid = isgraphics(hLines);
    legend(hLines(valid), {labelMap.(genotypes{1}), labelMap.(genotypes{2})}, ...
           'Location','southwest');
end

        end

        % Save one figure per condition
        safeCond = regexprep(cond, '[^a-zA-Z0-9]', '');
        outName  = fullfile(baseOutDir, sprintf('GroupPSD_%s.png', safeCond));
        saveas(gcf, outName);
        fprintf('Saved group PSD figure: %s\n', outName);
    end
end
