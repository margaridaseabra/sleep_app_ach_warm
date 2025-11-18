function OUT = run_bout_window_stats(PERHOUR, out_dir)
% RUN_BOUT_WINDOW_STATS
% -------------------------------------------------------------------------
% Compute bout-duration statistics in 3 time windows:
%   Window 1 : 0–3 h
%   Window 2 : 3–6 h (manipulation)
%   Window 3 : >6 h
%
% Using PERHOUR table (same as used for plotting), we build per-mouse
% mean bout duration per state×condition×genotype×window, then compute:
%
% 1) Repeated-measures ANOVA (Window × Genotype) per state×condition.
%    -> boutwin_anova_summary.csv
%
% 2) Within-genotype paired t-test: manipulation (W2) vs rest (mean of W1 & W3).
%    -> boutwin_manip_vs_rest_within_genotype.csv
%
% 3) Between-genotype unpaired t-tests: APP vs WT in each window (W1,W2,W3).
%    -> boutwin_genotype_diff_by_window.csv
%
% All tables are saved in out_dir and returned in OUT.
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

PH = PERHOUR;

needed = {'hour_idx','dur_s','bouts_per_h','state','condition','mouse','genotype'};
if ~all(ismember(needed, PH.Properties.VariableNames))
    error('PERHOUR table missing required columns for run_bout_window_stats.');
end

% ---- define windows consistent with the plots (0–3 h, 3–6 h, >6 h) -----
h = double(PH.hour_idx);
win_idx = zeros(size(h));           % 1,2,3
win_idx(h <= 2)          = 1;       % 0–3 h
win_idx(h >= 3 & h <= 5) = 2;       % 3–6 h
win_idx(h >= 6)          = 3;       % >6 h
PH.win = win_idx;

PH = PH(PH.win>=1 & PH.win<=3, :);
if isempty(PH)
    warning('run_bout_window_stats: no data in 0–3/3–6/>6 h windows.');
    OUT = struct();
    return;
end

% ---- normalized text columns -------------------------------------------
PH.state     = string(PH.state);
PH.condition = string(PH.condition);
PH.mouse     = string(PH.mouse);
PH.genotype  = string(PH.genotype);

geno = PH.genotype;
geno(geno ~= "WT") = "APP";   % collapse to WT vs APP
PH.geno_group = geno;

% ---- per-mouse mean bout duration per state×cond×geno×window ----------
[gid, st, cond, geno_g, mouse, win] = findgroups( ...
    PH.state, PH.condition, PH.geno_group, PH.mouse, PH.win);

tot_dur   = splitapply(@(x) sum(x,'omitnan'), double(PH.dur_s),      gid);
tot_bouts = splitapply(@(x) sum(x,'omitnan'), double(PH.bouts_per_h),gid);

mean_bout = tot_dur ./ max(tot_bouts,1);
mean_bout(tot_bouts==0) = NaN;

Tmouse = table(st, cond, geno_g, mouse, win, mean_bout, ...
    'VariableNames',{'state','condition','geno','mouse','win','mean_bout_dur_win_s'});

% ---- state / condition sets --------------------------------------------
state_pref = ["WK","NREM","REM"];
cond_pref  = ["baseline","ambtemp","drugs"];

states = state_pref(ismember(state_pref, unique(Tmouse.state)));
conds  = cond_pref(ismember(cond_pref, unique(Tmouse.condition)));

% ---- result rows (cell arrays) -----------------------------------------
anova_rows  = {};   % state, condition, p_genotype, p_window, p_interaction
manip_rows  = {};   % state, condition, genotype, n, mean_rest, mean_manip, diff, p_paired
geno_rows   = {};   % state, condition, window_label, n_WT, n_APP, mean_WT, mean_APP, diff, p_unpaired

for si = 1:numel(states)
    st_name = states(si);
    Tst = Tmouse(Tmouse.state == st_name, :);
    if isempty(Tst), continue; end

    conds_st = conds(ismember(conds, unique(Tst.condition)));
    if isempty(conds_st), continue; end

    for ci = 1:numel(conds_st)
        cnd = conds_st(ci);

        % ---- subtable for this state+condition -------------------------
        Tc = Tst(Tst.condition == cnd & ~isnan(Tst.mean_bout_dur_win_s), :);
        if height(Tc) < 3, continue; end

        % ---- build wide table manually: one row per mouse, cols W1,W2,W3
        umice = unique(Tc.mouse);
        nSub  = numel(umice);
        W1 = NaN(nSub,1);
        W2 = NaN(nSub,1);
        W3 = NaN(nSub,1);
        gsub = strings(nSub,1);

        for k = 1:nSub
            mk = umice(k);
            rowsM = Tc.mouse == mk;

            gk = unique(Tc.geno(rowsM));
            if numel(gk) ~= 1
                % something weird, skip this mouse
                continue;
            end
            gsub(k) = gk;

            % window 1..3
            for w = 1:3
                rowsMW = rowsM & Tc.win == w;
                if any(rowsMW)
                    val = Tc.mean_bout_dur_win_s(rowsMW);
                    mval = mean(val,'omitnan');
                    if w==1
                        W1(k) = mval;
                    elseif w==2
                        W2(k) = mval;
                    else
                        W3(k) = mval;
                    end
                end
            end
        end

        Tw = table(umice, gsub, W1, W2, W3, ...
            'VariableNames',{'mouse','geno','W1','W2','W3'});

        % drop rows with all-NaN windows
        allNaN = isnan(Tw.W1) & isnan(Tw.W2) & isnan(Tw.W3);
        Tw(allNaN,:) = [];
        if height(Tw) < 3, continue; end

        % ---- (1) RM ANOVA: Window × Genotype --------------------------
        pG = NaN; pW = NaN; pInt = NaN;

        measVars = {'W1','W2','W3'};
        % drop window columns that are all-NaN
        colAllNaN = all(isnan(Tw{:,measVars}),1);
        measVars  = measVars(~colAllNaN);

        if numel(measVars) >= 2
            % windows indices 1,2,3 for W1,W2,W3
            wnums = zeros(numel(measVars),1);
            for mv = 1:numel(measVars)
                wnums(mv) = str2double(measVars{mv}(2));  % 'W1'->1
            end
            WithinDesign = table(wnums,'VariableNames',{'Window'});

            formula = sprintf('%s-%s ~ geno', measVars{1}, measVars{end});

            try
                rm = fitrm(Tw, formula, 'WithinDesign', WithinDesign);

                % within-subject: Window and Window×geno
                rtbl = ranova(rm, 'WithinModel','Window');
                rn = rtbl.Properties.RowNames;
                rowW   = strcmp(rn,'Window');
                rowInt = strcmp(rn,'Window:geno');
                if any(rowW),   pW   = rtbl.pValue(rowW);   end
                if any(rowInt), pInt = rtbl.pValue(rowInt); end

                % between-subject: geno main effect
                bt = anova(rm);
                pG = NaN;
                if ismember('Term', bt.Properties.VariableNames)
                    term = string(bt.Term);
                    rowG = term=="geno";
                    if any(rowG) && ismember('pValue', bt.Properties.VariableNames)
                        pG = bt.pValue(rowG);
                    end
                else
                    rn_bt = bt.Properties.RowNames;
                    rowG = strcmp(rn_bt,'geno');
                    if any(rowG) && ismember('pValue', bt.Properties.VariableNames)
                        pG = bt.pValue(rowG);
                    end
                end
            catch
                % leave pG, pW, pInt as NaN if RM fails
            end
        end

        anova_rows(end+1,:) = {char(st_name), char(cnd), pG, pW, pInt}; %#ok<AGROW>

        % ---- (2) Within-genotype: manipulation vs rest -----------------
        % We need at least W2 and (W1 or W3). Here we use W1 & W3 both, like in plots.
        for gname = ["WT","APP"]
            maskG = Tw.geno == gname;
            if ~any(maskG), continue; end

            w1g = Tw.W1(maskG);
            w2g = Tw.W2(maskG);
            w3g = Tw.W3(maskG);

            rest  = nanmean([w1g, w3g], 2);   % mean of W1 and W3
            manip = w2g;

            valid = ~isnan(rest) & ~isnan(manip);
            if sum(valid) < 2, continue; end

            [~, p_paired] = ttest(manip(valid), rest(valid));

            mean_rest  = mean(rest(valid),'omitnan');
            mean_manip = mean(manip(valid),'omitnan');
            diff_val   = mean_manip - mean_rest;

            manip_rows(end+1,:) = { ...
                char(st_name), char(cnd), char(gname), ...
                sum(valid), mean_rest, mean_manip, diff_val, p_paired}; %#ok<AGROW>
        end

        % ---- (3) Between-genotype: APP vs WT per window ----------------
        for w = 1:3
            colname = sprintf('W%d', w);
            vals    = Tw.(colname);
            geno_tw = Tw.geno;

            WT_vals  = vals(geno_tw=="WT");
            APP_vals = vals(geno_tw=="APP");

            WT_vals  = WT_vals(~isnan(WT_vals));
            APP_vals = APP_vals(~isnan(APP_vals));

            if numel(WT_vals) < 2 || numel(APP_vals) < 2
                continue;
            end

            [~, p_unp] = ttest2(WT_vals, APP_vals);

            mean_WT  = mean(WT_vals,'omitnan');
            mean_APP = mean(APP_vals,'omitnan');
            diff_ga  = mean_APP - mean_WT;

            if w==1
                wlabel = '0–3 h';
            elseif w==2
                wlabel = '3–6 h (manip)';
            else
                wlabel = '>6 h';
            end

            geno_rows(end+1,:) = { ...
                char(st_name), char(cnd), wlabel, ...
                numel(WT_vals), numel(APP_vals), ...
                mean_WT, mean_APP, diff_ga, p_unp}; %#ok<AGROW>
        end
    end
end

% ---- Pack tables safely (handle empty cases) ---------------------------

% 1) ANOVA summary
if isempty(anova_rows)
    anova_tbl = table('Size',[0 5], ...
        'VariableTypes',{'string','string','double','double','double'}, ...
        'VariableNames',{'state','condition','p_genotype','p_window','p_interaction'});
else
    anova_tbl = cell2table(anova_rows, ...
        'VariableNames',{'state','condition','p_genotype','p_window','p_interaction'});
end

% 2) Manip vs rest (within genotype)
if isempty(manip_rows)
    manip_tbl = table('Size',[0 8], ...
        'VariableTypes',{'string','string','string','double','double','double','double','double'}, ...
        'VariableNames',{'state','condition','genotype','n', ...
                         'mean_rest','mean_manip','diff_manip_minus_rest','p_paired'});
else
    manip_tbl = cell2table(manip_rows, ...
        'VariableNames',{'state','condition','genotype','n', ...
                         'mean_rest','mean_manip','diff_manip_minus_rest','p_paired'});
end

% 3) Genotype diff by window
if isempty(geno_rows)
    geno_tbl = table('Size',[0 9], ...
        'VariableTypes',{'string','string','string', ...
                         'double','double','double','double','double','double'}, ...
        'VariableNames',{'state','condition','window', ...
                         'n_WT','n_APP','mean_WT','mean_APP', ...
                         'diff_APP_minus_WT','p_unpaired'});
else
    geno_tbl = cell2table(geno_rows, ...
        'VariableNames',{'state','condition','window', ...
                         'n_WT','n_APP','mean_WT','mean_APP', ...
                         'diff_APP_minus_WT','p_unpaired'});
end

% ---- Save to CSV -------------------------------------------------------
f1 = fullfile(out_dir, 'boutwin_anova_summary.csv');
f2 = fullfile(out_dir, 'boutwin_manip_vs_rest_within_genotype.csv');
f3 = fullfile(out_dir, 'boutwin_genotype_diff_by_window.csv');

writetable(anova_tbl, f1);
writetable(manip_tbl, f2);
writetable(geno_tbl, f3);

OUT = struct('anova',anova_tbl, ...
             'manip_vs_rest',manip_tbl, ...
             'geno_by_window',geno_tbl, ...
             'files',struct('anova_csv',f1, ...
                            'manip_csv',f2, ...
                            'geno_csv',f3));
end
