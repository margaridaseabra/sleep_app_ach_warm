function G2 = subset_group_by_cond(GROUP, cond_keep)
% G2 = subset_group_by_cond(GROUP, cond_keep)
%   cond_keep: string or cellstr, e.g. 'baseline' or {'baseline','ambtemp'}

    if ischar(cond_keep) || isstring(cond_keep)
        cond_keep = cellstr(cond_keep);
    end
    cond_keep = lower(string(cond_keep));

    sess_conds = lower(string({GROUP.sessions.cond}));
    keep_mask  = ismember(sess_conds, cond_keep);

    G2 = GROUP;
    G2.sessions = GROUP.sessions(keep_mask);
    G2.out      = GROUP.out(keep_mask);

    if isfield(GROUP, 'metrics')
        G2.metrics = GROUP.metrics(keep_mask);
    end

    if isfield(GROUP, 'metrics_tbl')
        G2.metrics_tbl = GROUP.metrics_tbl(keep_mask, :);
    else
        G2.metrics_tbl = struct2table(G2.metrics);
    end
end
