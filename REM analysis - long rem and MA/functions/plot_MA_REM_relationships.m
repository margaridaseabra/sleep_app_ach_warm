function plot_MA_REM_relationships(ALL_REM, out_dir)
% plot_MA_REM_relationships
% -------------------------------------------------------------
% Scatter plots:
%   - n_MA_pre vs REM dur, colored by geno

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end

ALL_REM.geno = string(ALL_REM.geno);
ALL_REM.cond = string(ALL_REM.cond);

conds = unique(ALL_REM.cond,'stable');
COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

for iC = 1:numel(conds)
    selC = ALL_REM.cond==conds(iC);
    R = ALL_REM(selC,:);
    if isempty(R); continue; end
    
    figure('Color','w'); hold on;
    
    for geno = ["WT","APP"]
        selG = R.geno==geno;
        if ~any(selG), continue; end
        
        if geno=="WT"; col = COL_WT; else; col = COL_APP; end
        
        scatter(R.n_MA_pre(selG), R.dur_s(selG), 25, col, 'filled', ...
                'MarkerFaceAlpha',0.6); hold on;
    end
    
    xlabel('Number of MAs in pre-REM window');
    ylabel('REM bout duration (s)');
    title(sprintf('MA vs REM duration – %s', conds(iC)));
    legend({'WT','APP'}, 'Location','best');
    box off;
    
    saveas(gcf, fullfile(out_dir, sprintf('MA_vs_REM_%s.png', conds(iC))));
end
