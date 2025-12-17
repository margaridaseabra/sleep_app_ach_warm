function OUT = export_baseline_REM_survival_forPrism(ALL_REM, out_dir, varargin)
% export_baseline_REM_survival_forPrism
% -------------------------------------------------------------------------
% Create a Prism-ready CSV to plot REM bout survival curves (WT vs APP).
%
% INPUT
%   ALL_REM : table with at least:
%       - ALL_REM.cond   : condition label (e.g. "baseline")
%       - ALL_REM.geno   : genotype ("WT" or "APP")
%       - ALL_REM.dur_s  : REM bout duration in seconds
%
%   out_dir : folder where the CSV will be saved
%
% NAME–VALUE OPTIONS
%   'baseline_cond' : which condition to use (default "baseline")
%
% OUTPUT (struct OUT)
%   OUT.success     : logical
%   OUT.nWT_bouts   : # of baseline WT bouts used
%   OUT.nAPP_bouts  : # of baseline APP bouts used
%   OUT.csv_file    : full path to CSV
%   OUT.table       : table written to CSV
%
% CSV FORMAT (Prism survival table)
%   WT_Time, WT_Censored, APP_Time, APP_Censored
%   12.3,    0,           10.1,     0
%   8.4,     0,           7.9,      0
%   ...      ...          ...       ...
%
% Notes:
%   - Each row is one bout (WT and APP columns can have NaN for padding).
%   - Censored flags are all 0 (no censoring; every bout ends).
% -------------------------------------------------------------------------

    p = inputParser;
    addParameter(p,'baseline_cond',"baseline",@(x)ischar(x) || isstring(x));
    parse(p, varargin{:});
    baseline_cond = string(p.Results.baseline_cond);

    if nargin < 2 || isempty(out_dir)
        out_dir = pwd;
    end
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    OUT = struct('success',false, ...
                 'nWT_bouts',0, ...
                 'nAPP_bouts',0, ...
                 'csv_file','', ...
                 'table',table());

    % ---- Basic checks ----
    if isempty(ALL_REM) || ~istable(ALL_REM)
        warning('ALL_REM is empty or not a table. Nothing to export.');
        return;
    end

    requiredVars = {'cond','geno','dur_s'};
    if ~all(ismember(requiredVars, ALL_REM.Properties.VariableNames))
        error('ALL_REM must contain variables: cond, geno, dur_s.');
    end

    % ---- Filter to baseline condition ----
    cond  = string(ALL_REM.cond);
    geno  = string(ALL_REM.geno);
    dur_s = ALL_REM.dur_s;

    mask_base = (cond == baseline_cond);
    if ~any(mask_base)
        warning('No rows in ALL_REM found for baseline condition "%s".', baseline_cond);
        return;
    end

    dur_base  = dur_s(mask_base);
    geno_base = geno(mask_base);

    % ---- Split by genotype, clean durations ----
    dWT  = dur_base(geno_base == "WT");
    dAPP = dur_base(geno_base == "APP");

    dWT  = dWT(~isnan(dWT) & dWT > 0);
    dAPP = dAPP(~isnan(dAPP) & dAPP > 0);

    nWT  = numel(dWT);
    nAPP = numel(dAPP);

    if nWT == 0 && nAPP == 0
        warning('No positive REM bout durations found for WT or APP in baseline "%s".', baseline_cond);
        return;
    end

    % ---- Build Prism-style table ----
    maxN = max(nWT, nAPP);

    WT_Time      = nan(maxN,1);
    WT_Censored  = nan(maxN,1);
    APP_Time     = nan(maxN,1);
    APP_Censored = nan(maxN,1);

    if nWT > 0
        WT_Time(1:nWT)     = dWT;
        WT_Censored(1:nWT) = 0;   % 0 = event, not censored
    end
    if nAPP > 0
        APP_Time(1:nAPP)     = dAPP;
        APP_Censored(1:nAPP) = 0; % 0 = event, not censored
    end

    Tprism = table(WT_Time, WT_Censored, APP_Time, APP_Censored);

    % ---- Write CSV ----
    % Make a safe filename based on the condition
    cond_tag  = matlab.lang.makeValidName(char(baseline_cond),'ReplacementStyle','delete');
    if isempty(cond_tag)
        cond_tag = 'baseline';
    end

    csv_name = sprintf('REM_survival_%s_forPrism.csv', cond_tag);
    csv_file = fullfile(out_dir, csv_name);

    try
        writetable(Tprism, csv_file);
        fprintf('📄 Prism survival CSV saved to: %s\n', csv_file);
    catch ME
        warning('export_baseline_REM_survival_forPrism:WriteFailed', ...
                'Could not write CSV (%s): %s', csv_file, ME.message);
        return;
    end

    % ---- Pack OUT ----
    OUT.success    = true;
    OUT.nWT_bouts  = nWT;
    OUT.nAPP_bouts = nAPP;
    OUT.csv_file   = csv_file;
    OUT.table      = Tprism;
end
