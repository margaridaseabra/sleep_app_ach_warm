function GROUP_TRANS = run_ach_peri_from_GROUP(groupMatFile)
% RUN_ACH_PERI_FROM_GROUP
%   Use existing ACh_GROUP_all.mat to run peri-onset ACh analysis
%   for ALL sessions (Wake, NREM, REM onsets).
%
%   GROUP_TRANS is saved to ACh_periTransitions_all.mat in same folder.

if nargin < 1 || isempty(groupMatFile)
    [fname, fpath] = uigetfile('*.mat', 'Select ACh_GROUP_all.mat');
    if isequal(fname,0)
        error('No file selected.');
    end
    groupMatFile = fullfile(fpath, fname);
end

S = load(groupMatFile);
if ~isfield(S,'GROUP')
    error('Selected .mat file does not contain GROUP struct.');
end
GROUP = S.GROUP;

% Make sure codes match what you used in scoring / ach_analysis
CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);

nSess = numel(GROUP.sessions);
GROUP_TRANS.sessions = GROUP.sessions;
GROUP_TRANS.out      = cell(nSess,1);
GROUP_TRANS.codes    = CODES;

out_root = fileparts(groupMatFile);

fprintf('\nPeri-onset ACh analysis for %d sessions.\n', nSess);

for k = 1:nSess
    s = GROUP.sessions(k);
    fprintf('\n=== %d/%d: %s | %s | %s ===\n', ...
            k, nSess, s.mouse, s.geno, s.cond);

    label     = sprintf('%s_%s_%s', s.mouse, s.geno, s.cond);
    out_pref  = label;

    OUT_T = ach_peri_transitions_single( ...
        s.mat, s.csv, CODES, ...
        't_pre', 30, ...
        't_post', 60, ...
        'session_label', label, ...
        'out_prefix', out_pref, ...
        'out_dir', out_root, ...
        'zscore_rows', false, ...
        'save_fig', true);

    GROUP_TRANS.out{k} = OUT_T;
end

save(fullfile(out_root, 'ACh_periTransitions_all.mat'), 'GROUP_TRANS');
fprintf('\nSaved ACh_periTransitions_all.mat in %s\n', out_root);
end
