input_dir = '/Users/margaridaseabra/24.11scores';


OUT_REM = run_REM_length_structure_analysis(input_dir, ...
    'pattern', '*_scores_1Hz.csv', ...
    'codes', struct('WK',0,'NREM',1,'REM',2,'MA',15), ...
    'quantiles', [33 66], ...       % defines short / medium / long
    'cluster_gap_sec', 600, ...     % 10 min max gap inside a REM cluster
    'out_dir', fullfile(input_dir,'rem_analysis'));


