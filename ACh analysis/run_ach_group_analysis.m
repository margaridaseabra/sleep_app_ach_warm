function GROUP = run_ach_group_analysis()

% ---- scoring codes (adjust if needed) ----
CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);

% --------- DEFINE ALL SESSIONS HERE -----------------
% One row per .mat / condition / mouse
% Fill with your real paths:
i = 0;  SESS = struct([]);

% WT BASELINE
i=i+1; SESS(i) = struct( ...
    'mouse','M1', ...
    'geno','WT', ...
    'cond','baseline', ...
    'mat','/path/to/M1_base.mat', ...
    'csv','/path/to/M1_base_scores_1Hz.csv');

% WT AMBTEMP
i=i+1; SESS(i) = struct( ...
    'mouse','M1', ...
    'geno','WT', ...
    'cond','ambtemp', ...
    'mat','/path/to/M1_ambtemp.mat', ...
    'csv','/path/to/M1_ambtemp_scores_1Hz.csv');

% WT DRUG
i=i+1; SESS(i) = struct( ...
    'mouse','M1', ...
    'geno','WT', ...
    'cond','drug', ...
    'mat','/path/to/M1_drug.mat', ...
    'csv','/path/to/M1_drug_scores_1Hz.csv');

% ---- APP/PS1 sessions for M2, etc...
% i=i+1; SESS(i) = struct('mouse','M2','geno','APPPS1',...
%                         'cond','baseline',...
%                         'mat','...', 'csv','...');

nSess = numel(SESS);

% ------------- RUN ach_analysis FOR EACH SESSION ---------------
GROUP.sessions = SESS;
GROUP.out      = cell(nSess,1);

for k = 1:nSess
    s = SESS(k);
    fprintf('\n=== Processing %s %s (%s) ===\n', s.mouse, s.geno, s.cond);

    OUT = ach_analysis( ...
        s.mat, s.csv, ...
        'codes', CODES, ...
        'mouse_id', s.mouse, ...
        'session', s.cond, ...
        'out_prefix', sprintf('%s_%s_%s', s.mouse, s.geno, s.cond), ...
        'verbose', false);

    GROUP.out{k} = OUT;

    % Store a light summary row directly in GROUP.metrics
    GROUP.metrics(k).mouse      = s.mouse;
    GROUP.metrics(k).geno       = s.geno;
    GROUP.metrics(k).cond       = s.cond;

    % NREM PSD summary
    if ~isempty(OUT.psd.NREM.f)
        GROUP.metrics(k).NREM_power    = OUT.psd.NREM.band_power;
        GROUP.metrics(k).NREM_peakHz   = OUT.psd.NREM.peak_freq;
        GROUP.metrics(k).NREM_peakAmp  = OUT.psd.NREM.peak_amp;
    else
        GROUP.metrics(k).NREM_power    = NaN;
        GROUP.metrics(k).NREM_peakHz   = NaN;
        GROUP.metrics(k).NREM_peakAmp  = NaN;
    end

    % State slopes
    ss = OUT.state_slopes;
    GROUP.metrics(k).slope_Wake = getfield_or_nan(ss,'Wake');
    GROUP.metrics(k).slope_NREM = getfield_or_nan(ss,'NREM');
    GROUP.metrics(k).slope_REM  = getfield_or_nan(ss,'REM');

    % Wake_onset transition slopes & peaks (if present)
    idxW = find(strcmp({OUT.transitions.name},'Wake_onset'));
    if ~isempty(idxW)
        Tw = OUT.transitions(idxW);
        GROUP.metrics(k).WakeOn_peak_mean  = mean(Tw.peaks,'omitnan');
        GROUP.metrics(k).WakeOn_slope_mean = mean(Tw.slopes,'omitnan');
    else
        GROUP.metrics(k).WakeOn_peak_mean  = NaN;
        GROUP.metrics(k).WakeOn_slope_mean = NaN;
    end
end

% Turn the metrics struct into a table and save
GROUP.metrics_tbl = struct2table(GROUP.metrics);
writetable(GROUP.metrics_tbl,'ACh_group_metrics.csv');

% ------------- PLOTS -----------------
group_plot_psd_NREM(GROUP);
group_plot_wake_onsets(GROUP);

save('ACh_GROUP_all.mat','GROUP');
end

% small helper:
function m = getfield_or_nan(ss, field)
if isfield(ss,field)
    m = ss.(field).mean;
else
    m = NaN;
end
end
