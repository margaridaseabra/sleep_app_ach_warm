function group_plot_all_transitions(GROUP, out_dir, show_mouse_ids)
% Make paper-style transition figures with optional mouse ID labels
%
% INPUTS:
%   GROUP          - structure with sessions and transitions
%   out_dir        - directory to save figures (optional)
%   show_mouse_ids - logical, if true show mouse IDs on scatter plots (default: false)

if nargin < 2
    out_dir = [];
end
if nargin < 3 || isempty(show_mouse_ids)
    show_mouse_ids = false;
end

trans_list  = {'Wake_onset','NREM_onset','REM_onset'};
title_list  = {'WAKE ONSETS','NREM ONSETS','REM ONSETS'};

% === DIAGNOSTIC: Check what transitions are available ===
fprintf('\n=== Transition Data Availability ===\n');
sessions = GROUP.sessions;
for k = 1:numel(sessions)
    fprintf('Session %d (%s, %s, %s): ', k, sessions(k).mouse, sessions(k).geno, sessions(k).cond);
    OUT = GROUP.out{k};
    if isfield(OUT, 'transitions') && ~isempty(OUT.transitions)
        trans_names = {OUT.transitions.name};
        fprintf('%s\n', strjoin(trans_names, ', '));
    else
        fprintf('NO TRANSITIONS\n');
    end
end
fprintf('====================================\n\n');

for i = 1:numel(trans_list)
    fig = plot_one_transition(GROUP, trans_list{i}, title_list{i}, show_mouse_ids);

    if ~isempty(out_dir)
        if ~exist(out_dir,'dir'); mkdir(out_dir); end
        fname = fullfile(out_dir, sprintf('Transitions_%s.png', trans_list{i}));
        saveas(fig, fname);
    end
end
end

% ======================================================================
function fig = plot_one_transition(GROUP, trans_name, big_title, show_mouse_ids)

sessions = GROUP.sessions;
conds    = unique(string({sessions.cond}),'stable');
nCond    = numel(conds);

% figure out genotype labels and order (WT first if present)
geno_vals = unique(string({sessions.geno}),'stable');
if any(geno_vals=="WT")
    genoWT  = "WT";
    genoMut = geno_vals(geno_vals~="WT");
else
    genoWT  = geno_vals(1);
    genoMut = geno_vals(2:end);
end
if isempty(genoMut)
    error('Need at least two genotypes (WT + APP).');
end
genoMut = genoMut(1);   % assume one mutant group (APP/PS1)

% colours (WT grey, APP cornflower blue)
COL_WT  = [0.6 0.6 0.6];           % grey
COL_APP = [0.392 0.584 0.929];     % cornflower blue

% === DIAGNOSTIC: Check conditions and genotypes ===
fprintf('\n=== %s - Conditions and Genotypes ===\n', trans_name);
fprintf('Conditions found: %s\n', strjoin(conds, ', '));
fprintf('WT genotype: %s\n', genoWT);
fprintf('Mutant genotype: %s\n', genoMut);

% Check distribution
for c = 1:numel(conds)
    n_wt = sum(string({sessions.cond}) == conds(c) & string({sessions.geno}) == genoWT);
    n_app = sum(string({sessions.cond}) == conds(c) & string({sessions.geno}) == genoMut);
    fprintf('  %s: WT n=%d, %s n=%d\n', conds(c), n_wt, genoMut, n_app);
end
fprintf('=====================================\n\n');

fig = figure('Name',sprintf('%s – %s',big_title,trans_name),'Color','w', ...
             'Position',[100 100 1200 300*nCond]);

for c = 1:nCond
    cond = conds(c);

    WT_traces   = [];
    APP_traces  = [];
    WT_peaks    = [];
    APP_peaks   = [];
    WT_mice     = {};
    APP_mice    = {};
    t_ref       = [];

    % ---------- collect per-session means for this condition ----------
    for k = 1:numel(sessions)
        if string(sessions(k).cond) ~= cond
            continue;
        end

        OUT = GROUP.out{k};

        % find this transition in OUT.transitions
        idxT = [];
        if isfield(OUT,'transitions') && ~isempty(OUT.transitions)
            idxT = find(strcmp({OUT.transitions.name}, trans_name));
        end
        if isempty(idxT)
            continue;
        end

        Ttr = OUT.transitions(idxT);
        
        % Check if we have valid mean trace
        if ~isfield(Ttr, 'mean') || isempty(Ttr.mean) || all(isnan(Ttr.mean))
            continue;
        end
        
        % Get time axis - check different possible field names
        if isfield(Ttr, 't_axis')
            t_this = Ttr.t_axis(:);
        elseif isfield(Ttr, 'time')
            t_this = Ttr.time(:);
        elseif isfield(Ttr, 't')
            t_this = Ttr.t(:);
        else
            % Create time axis from window info if available
            if isfield(Ttr, 'win_sec')
                n_pts = numel(Ttr.mean);
                t_this = linspace(Ttr.win_sec(1), Ttr.win_sec(2), n_pts)';
            else
                % Default: assume symmetric window
                n_pts = numel(Ttr.mean);
                t_this = linspace(-100, 100, n_pts)';
            end
        end

        % reference time axis
        if isempty(t_ref)
            t_ref = t_this;
        elseif numel(t_this) ~= numel(t_ref) || any(abs(t_this - t_ref) > 1e-6)
            % interpolate to reference
            Ttr.mean = interp1(t_this, Ttr.mean, t_ref, 'linear', NaN);
        end

        g = string(sessions(k).geno);
        
        % --------------------------
        %     PEAK COMPUTATION (same as temperature version)
        % --------------------------
        if isfield(Ttr, 'peak_mean') && ~isnan(Ttr.peak_mean)
            peak_val = Ttr.peak_mean;
        elseif isfield(Ttr, 'peak') && ~isnan(Ttr.peak)
            peak_val = Ttr.peak;
        else
            % Compute peak from mean trace
            peak_val = max(Ttr.mean(:), [], 'omitnan');
        end
        
        % Add to appropriate group
        if g == genoWT
            WT_traces = [WT_traces Ttr.mean(:)]; %#ok<AGROW>
            if ~isnan(peak_val)
                WT_peaks = [WT_peaks; peak_val]; %#ok<AGROW>
                WT_mice{end+1} = sessions(k).mouse; %#ok<AGROW>
            end
        elseif g == genoMut
            APP_traces = [APP_traces Ttr.mean(:)]; %#ok<AGROW>
            if ~isnan(peak_val)
                APP_peaks = [APP_peaks; peak_val]; %#ok<AGROW>
                APP_mice{end+1} = sessions(k).mouse; %#ok<AGROW>
            end
        end
    end

    fprintf('%s - %s: WT=%d traces, APP=%d traces, WT peaks=%d, APP peaks=%d\n', ...
            char(cond), trans_name, size(WT_traces,2), size(APP_traces,2), ...
            numel(WT_peaks), numel(APP_peaks));

    % ========== LEFT: mean trace ± SEM ==================
    subplot(nCond, 2, 2*c-1); hold on;
    
    if ~isempty(WT_traces) && ~all(isnan(WT_traces(:)))
        m_wt  = mean(WT_traces, 2, 'omitnan');
        se_wt = std(WT_traces, 0, 2, 'omitnan') / sqrt(size(WT_traces,2));
        fill_between(t_ref, m_wt - se_wt, m_wt + se_wt, COL_WT, 0.2);
        plot(t_ref, m_wt, 'Color', COL_WT, 'LineWidth', 1.5);
    end
    
    if ~isempty(APP_traces) && ~all(isnan(APP_traces(:)))
        m_app  = mean(APP_traces, 2, 'omitnan');
        se_app = std(APP_traces, 0, 2, 'omitnan') / sqrt(size(APP_traces,2));
        fill_between(t_ref, m_app - se_app, m_app + se_app, COL_APP, 0.2);
        plot(t_ref, m_app, 'Color', COL_APP, 'LineWidth', 1.5);
    end
    
    plot([0 0], ylim, 'k--', 'LineWidth', 1);
    xlabel('Time from transition (s)');
    ylabel('ACh (ΔF/F)');
    title(sprintf('%s', cond));
    box off;
    if c == 1 && (~isempty(WT_traces) || ~isempty(APP_traces))
        legend({char(genoWT), char(genoMut)}, 'Box', 'off', 'Location', 'best');
    end

    % ========== RIGHT: peak bars + individual points ==================
    subplot(nCond, 2, 2*c); hold on;

    bar(1, mean(WT_peaks,  'omitnan'), 0.6, 'FaceColor', COL_WT);
    bar(2, mean(APP_peaks, 'omitnan'), 0.6, 'FaceColor', COL_APP);

    % scatter individual sessions
    jitter = 0.06;
    if ~isempty(WT_peaks)
        xpos = 1 + (rand(size(WT_peaks))-0.5)*2*jitter;
        scatter(xpos, WT_peaks, 25, COL_WT,'filled');
        
        % Add mouse IDs if requested
        if show_mouse_ids
            for i = 1:numel(WT_peaks)
                if i <= numel(WT_mice)
                    text(xpos(i), WT_peaks(i), sprintf(' %s', WT_mice{i}), ...
                         'FontSize', 7, 'HorizontalAlignment', 'left');
                end
            end
        end
    end
    if ~isempty(APP_peaks)
        xpos = 2 + (rand(size(APP_peaks))-0.5)*2*jitter;
        scatter(xpos, APP_peaks, 25, COL_APP,'filled');
        
        % Add mouse IDs if requested
        if show_mouse_ids
            for i = 1:numel(APP_peaks)
                if i <= numel(APP_mice)
                    text(xpos(i), APP_peaks(i), sprintf(' %s', APP_mice{i}), ...
                         'FontSize', 7, 'HorizontalAlignment', 'left');
                end
            end
        end
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{char(genoWT), char(genoMut)});
    ylabel('Peak ACh (ΔF/F)');
    title(sprintf('%s – peak response', cond));
    box off;

    % t-test with p-value display
    if ~isempty(WT_peaks) && ~isempty(APP_peaks)
        [~, p_val] = ttest2(WT_peaks, APP_peaks);
        star = p_to_star(p_val);
        
        % Display statistics
        ymax = max([WT_peaks(:); APP_peaks(:)]);
        ymin = min([WT_peaks(:); APP_peaks(:)]);
        yrange = ymax - ymin;
        
        if ~isempty(star)
            % Show significance stars
            ypos = ymax + yrange * 0.1;
            plot([1 2], [ypos ypos], 'k-', 'LineWidth', 1);
            text(1.5, ypos, star, 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'bold');
        end
        
        % Always show p-value below the plot
        if p_val < 0.001
            p_text = 'p < 0.001';
        else
            p_text = sprintf('p = %.3f', p_val);
        end
        text(1.5, ymin - yrange * 0.15, p_text, ...
             'HorizontalAlignment', 'center', 'FontSize', 9);
        
        % Also show means ± SEM
        fprintf('  %s %s - WT: %.3f±%.3f, APP: %.3f±%.3f, p=%.4f %s\n', ...
                char(cond), trans_name, ...
                mean(WT_peaks,'omitnan'), std(WT_peaks,'omitnan')/sqrt(numel(WT_peaks)), ...
                mean(APP_peaks,'omitnan'), std(APP_peaks,'omitnan')/sqrt(numel(APP_peaks)), ...
                p_val, star);
    elseif isempty(WT_peaks) || isempty(APP_peaks)
        % Note when comparison cannot be done
        text(1.5, mean(ylim), 'Insufficient data', ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.5 0.5 0.5]);
    end
end

sgtitle(big_title);

end

% ----------------------------------------------------------------------
function fill_between(x,y1,y2,col,alphaVal)
% small helper for shaded SEM
x   = x(:);
y1  = y1(:);
y2  = y2(:);
fill([x; flipud(x)], [y1; flipud(y2)], col, ...
     'EdgeColor','none','FaceAlpha',alphaVal);
end

% ----------------------------------------------------------------------
function star = p_to_star(p)
% Convert p-value to significance stars
if p < 0.001
    star = '***';
elseif p < 0.01
    star = '**';
elseif p < 0.05
    star = '*';
else
    star = '';
end
end