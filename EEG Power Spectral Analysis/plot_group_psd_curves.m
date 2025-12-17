function plot_group_psd_curves(baseOutDir)
% plot_group_psd_curves(baseOutDir)
%
% baseOutDir: root folder with structure
%   baseOutDir/genotype/condition/mouseID/PSD_data.mat
%
% PART 1 (unchanged logic, slightly renamed outputs):
%   For each CONDITION, creates a figure with 3 subplots:
%       Wake, NREM, REM
%   Each subplot shows mean PSD (in dB) for WT vs APP with SEM shading.
%
% PART 2 (NEW):
%   For each GENOTYPE, creates a figure with 3 subplots:
%       Wake, NREM, REM
%   Each subplot shows mean PSD (in dB) for *conditions* (e.g. baseline vs ambtemp)
%   within that genotype, with SEM shading.

    if nargin < 1 || isempty(baseOutDir)
        baseOutDir = 'EEG_PSD_AllMice_AmbTemp';
    end

    % ------------ Genotypes & colours ------------
    genotypes = {'WT','APP'};   % change here if your folders are different

    % Labels for legend (by genotype):
    labelMap = struct('WT','WT', 'APP','APP/PS1');

    % Colours for genotype comparison (WT grey, APP cornflower blue)
    colMap.WT  = [0.6 0.6 0.6];
    colMap.APP = [0.392 0.584 0.929];

    % States to plot
    states = {'Wake','NREM','REM'};

    % Common frequency axis (Hz) for all mice
    F_common = (0:0.25:100)';   % 0–100 Hz in 0.25 Hz steps

    % ------------ Collect all conditions ------------
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

    %% ============================================================
    %  PART 1: Compare GENOTYPES within each CONDITION (original view)
    % ==============================================================
    for ci = 1:numel(conditions)
        cond = conditions{ci};

        figure('Color','w','Name',cond);
        nStates = numel(states);

        for si = 1:nStates
            st = states{si};
            subplot(1, nStates, si); hold on;

            % Store line handles for legend
            hLines = gobjects(1, numel(genotypes));

            % --------- For each genotype, collect PSDs ----------
            for gi = 1:numel(genotypes)
                geno = genotypes{gi};
                gDir = fullfile(baseOutDir, geno, cond);
                if ~exist(gDir,'dir')
                    warning('Missing folder: %s', gDir);
                    continue;
                end

                mList = dir(gDir);
                mList = mList([mList.isdir] & ~ismember({mList.name},{'.','..'}));

                allPSD_dB = [];  % each column = one mouse

                for mi = 1:numel(mList)
                    sesDir  = fullfile(gDir, mList(mi).name);
                    matFile = fullfile(sesDir, 'PSD_data.mat');
                    if ~exist(matFile,'file')
                        warning('No PSD_data.mat in %s', sesDir);
                        continue;
                    end

                    S = load(matFile);
                    if ~isfield(S,'result') || ~isfield(S.result, st)
                        continue;
                    end

                    F   = S.result.(st).F;
                    PSD = S.result.(st).PSD;      % linear units

                    % Restrict to 0–100 Hz
                    idxBand = F >= 0 & F <= 100;
                    Fseg    = F(idxBand);
                    PSDseg  = PSD(idxBand);

                    % Convert to dB
                    PSDdB = 10*log10(PSDseg);

                    % Interpolate onto common frequency axis
                    PSDdB_common = interp1(Fseg, PSDdB, F_common, 'linear', 'extrap');

                    allPSD_dB(:, end+1) = PSDdB_common(:); %#ok<AGROW>
                end

                if isempty(allPSD_dB)
                    warning('No PSD data for %s, %s, %s', cond, st, geno);
                    continue;
                end

                % Mean and SEM across mice
                mPSD  = mean(allPSD_dB, 2, 'omitnan');
                nSubs = sum(~isnan(allPSD_dB(1,:)));
                sPSD  = std(allPSD_dB, 0, 2, 'omitnan') ./ max(nSubs,1).^0.5;

                c = colMap.(geno);

                x      = F_common;
                upper  = mPSD + sPSD;
                lower  = mPSD - sPSD;

                % Shaded SEM (very light)
                fill([x; flipud(x)], [upper; flipud(lower)], c, ...
                     'FaceAlpha', 0.20, 'EdgeColor', 'none');

                % Mean line (slightly darker + thicker)
                hLines(gi) = plot(x, mPSD, 'Color', c*0.9, 'LineWidth', 2.0);
            end

            % --------- Axis cosmetics ----------
            xlim([0 100]);
            xlabel('Frequency (Hz)');
            if si == 1
                ylabel('Power (dB)');
            end
            title(st);
            grid on;

            % Legend on the last subplot (REM)
            if si == nStates
                valid = isgraphics(hLines);
                legLabels = {};
                for gi = find(valid)
                    gname = genotypes{gi};
                    if isfield(labelMap, gname)
                        legLabels{end+1} = labelMap.(gname); %#ok<AGROW>
                    else
                        legLabels{end+1} = gname; %#ok<AGROW>
                    end
                end
                if ~isempty(legLabels)
                    legend(hLines(valid), legLabels, 'Location','southwest');
                end
            end
        end

        % Big title with condition name
        sgtitle(sprintf('Condition: %s (WT vs APP)', cond), 'FontWeight','bold');

        % Save figure for this condition
        safeCond = regexprep(cond, '[^a-zA-Z0-9]', '');
        outName  = fullfile(baseOutDir, sprintf('GroupPSD_%s_byGenotype.png', safeCond));
        saveas(gcf, outName);
        fprintf('Saved group PSD figure (by genotype): %s\n', outName);
    end

    %% ============================================================
    %  PART 2: Compare CONDITIONS within each GENOTYPE (new view)
    % ============================================================
    % Colour map for conditions (e.g. baseline vs ambtemp)
    condColors = lines(numel(conditions));  % one colour per condition

    for gi = 1:numel(genotypes)
        geno = genotypes{gi};

        figure('Color','w','Name',geno);
        nStates = numel(states);

        for si = 1:nStates
            st = states{si};
            subplot(1, nStates, si); hold on;

            % Store line handles for legend (one per condition)
            hLines = gobjects(1, numel(conditions));

            % --------- For each condition, collect PSDs ----------
            for ci = 1:numel(conditions)
                cond = conditions{ci};
                gDir = fullfile(baseOutDir, geno, cond);
                if ~exist(gDir,'dir')
                    warning('Missing folder: %s', gDir);
                    continue;
                end

                mList = dir(gDir);
                mList = mList([mList.isdir] & ~ismember({mList.name},{'.','..'}));

                allPSD_dB = [];  % each column = one mouse

                for mi = 1:numel(mList)
                    sesDir  = fullfile(gDir, mList(mi).name);
                    matFile = fullfile(sesDir, 'PSD_data.mat');
                    if ~exist(matFile,'file')
                        warning('No PSD_data.mat in %s', sesDir);
                        continue;
                    end

                    S = load(matFile);
                    if ~isfield(S,'result') || ~isfield(S.result, st)
                        continue;
                    end

                    F   = S.result.(st).F;
                    PSD = S.result.(st).PSD;      % linear units

                    % Restrict to 0–100 Hz
                    idxBand = F >= 0 & F <= 100;
                    Fseg    = F(idxBand);
                    PSDseg  = PSD(idxBand);

                    % Convert to dB
                    PSDdB = 10*log10(PSDseg);

                    % Interpolate onto common frequency axis
                    PSDdB_common = interp1(Fseg, PSDdB, F_common, 'linear', 'extrap');

                    allPSD_dB(:, end+1) = PSDdB_common(:); %#ok<AGROW>
                end

                if isempty(allPSD_dB)
                    warning('No PSD data for genotype %s, condition %s, state %s', geno, cond, st);
                    continue;
                end

                % Mean and SEM across mice
                mPSD  = mean(allPSD_dB, 2, 'omitnan');
                nSubs = sum(~isnan(allPSD_dB(1,:)));
                sPSD  = std(allPSD_dB, 0, 2, 'omitnan') ./ max(nSubs,1).^0.5;

                x      = F_common;
                upper  = mPSD + sPSD;
                lower  = mPSD - sPSD;

                c = condColors(ci,:);

                % Shaded SEM
                fill([x; flipud(x)], [upper; flipud(lower)], c, ...
                     'FaceAlpha', 0.20, 'EdgeColor', 'none');

                % Mean line
                hLines(ci) = plot(x, mPSD, 'Color', c*0.9, 'LineWidth', 2.0);
            end

            % --------- Axis cosmetics ----------
            xlim([0 100]);
            xlabel('Frequency (Hz)');
            if si == 1
                ylabel('Power (dB)');
            end
            title(st);
            grid on;

            % Legend on last subplot (REM)
            if si == nStates
                valid = isgraphics(hLines);
                if any(valid)
                    % Use the condition names directly as legend labels
                    legend(hLines(valid), conditions(valid), 'Location','southwest');
                end
            end
        end

        % Big title with genotype name
        sgtitle(sprintf('Genotype: %s (conditions compared)', geno), 'FontWeight','bold');

        % Save figure for this genotype
        safeGeno = regexprep(geno, '[^a-zA-Z0-9]', '');
        outName  = fullfile(baseOutDir, sprintf('GroupPSD_%s_byCondition.png', safeGeno));
        saveas(gcf, outName);
        fprintf('Saved group PSD figure (by condition): %s\n', outName);
    end
end
