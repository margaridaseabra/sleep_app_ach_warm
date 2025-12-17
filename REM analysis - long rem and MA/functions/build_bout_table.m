function BOUTS = build_bout_table(states, epochs_t, mouse_id, geno, cond)
% build_bout_table
% -------------------------------------------------------------
% From per-epoch state labels → bout table (W/NREM/REM).
%
% Inputs:
%   states   : [nEpochs x 1] (string, char, or numeric)
%   epochs_t : [nEpochs x 1] epoch start times (s)
%   mouse_id : scalar string/char
%   geno     : scalar string/char
%   cond     : scalar string/char

states = string(states(:));
epochs_t = epochs_t(:);
mouse_id = string(mouse_id);
geno     = string(geno);
cond     = string(cond);

nEpochs = numel(states);
if numel(epochs_t) ~= nEpochs
    error('states and epochs_t must have same length.');
end

% assume constant epoch length
if nEpochs > 1
    epoch_len_s = epochs_t(2) - epochs_t(1);
else
    epoch_len_s = 4;  % default
end

% find boundaries where state changes
is_new_bout = [true; states(2:end) ~= states(1:end-1)];
bout_start_idx = find(is_new_bout);

% bout end index (previous idx before next start, or last epoch)
bout_end_idx = [bout_start_idx(2:end)-1; nEpochs];

nBouts = numel(bout_start_idx);

mouse_col = strings(nBouts,1);
geno_col  = strings(nBouts,1);
cond_col  = strings(nBouts,1);
state_col = strings(nBouts,1);
t_start   = nan(nBouts,1);
t_end     = nan(nBouts,1);
dur_s     = nan(nBouts,1);

for iB = 1:nBouts
    iBeg = bout_start_idx(iB);
    iEnd = bout_end_idx(iB);
    
    mouse_col(iB) = mouse_id;
    geno_col(iB)  = geno;
    cond_col(iB)  = cond;
    state_col(iB) = states(iBeg);
    
    t_start(iB) = epochs_t(iBeg);
    t_end(iB)   = epochs_t(iEnd) + epoch_len_s;
    dur_s(iB)   = t_end(iB) - t_start(iB);
end

BOUTS = table;
BOUTS.mouse   = mouse_col;
BOUTS.geno    = geno_col;
BOUTS.cond    = cond_col;
BOUTS.state   = state_col;
BOUTS.t_start = t_start;
BOUTS.t_end   = t_end;
BOUTS.dur_s   = dur_s;
