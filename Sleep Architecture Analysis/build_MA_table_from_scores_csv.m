function MA = build_MA_table_from_scores_csv(scores_csv, mouse_id, genotype, condition, varargin)
% build_MA_table_from_scores_csv
% -------------------------------------------------------------------------
% Build a per-episode micro-arousal (MA) table from a 1 Hz scored CSV.
%
% For ONE recording:
%   - Detect contiguous runs of MA code (e.g. 15)
%   - For each MA bout, compute:
%       * onset time (s from recording start)
%       * offset time (s)
%       * state_before (NREM / REM / OTHER), based on the last non-MA
%         sample immediately before MA onset.
%
% Returned table MA has columns:
%   mouse        : string
%   genotype     : string
%   condition    : string
%   ma_on_s      : MA onset time (s)
%   ma_off_s     : MA offset time (s)
%   state_before : "NREM","REM","WK","OTHER","UNK"
%
% OPTIONS (name–value):
%   'codes'     : struct with fields WK,NREM,REM,MA (default: 1/4/9/15)
%   'epoch_sec' : sampling interval in seconds (default: 1)
%   'stateCol'  : name of column in CSV with numeric sleep codes
%                 (default: 'state' – change to 'score' or whatever you use)
%
% Usage example:
%   MA1 = build_MA_table_from_scores_csv('mouse1_baseline_1Hz.csv', ...
%           'mouse1','WT','baseline', ...
%           'codes', struct('WK',1,'NREM',4,'REM',9,'MA',15), ...
%           'stateCol','score');
% -------------------------------------------------------------------------

% --------- Parse options ----------
p = inputParser;
addParameter(p,'codes',struct('WK',1,'NREM',4,'REM',9,'MA',15),@isstruct);
addParameter(p,'epoch_sec',1,@(x)isscalar(x) && x>0);
addParameter(p,'stateCol','state',@(s)ischar(s) || isstring(s));
parse(p, varargin{:});

codes     = p.Results.codes;
epoch_sec = p.Results.epoch_sec;
stateCol  = string(p.Results.stateCol);

% --------- Read scores CSV ----------
T = readtable(scores_csv);

if ~ismember(stateCol, string(T.Properties.VariableNames))
    error('build_MA_table_from_scores_csv: column "%s" not found in %s. Adapt stateCol.', ...
        stateCol, scores_csv);
end

codes_vec = T.(stateCol);   % assume numeric codes (e.g. 0,1,4,9,15, ...)

if ~isnumeric(codes_vec)
    error('build_MA_table_from_scores_csv: state column "%s" is not numeric. Adapt code accordingly.', stateCol);
end

codes_vec = codes_vec(:);   % column vector

% --------- Detect MA bouts ----------
isMA = (codes_vec == codes.MA);

% add padding to catch rises/falls
d = diff([false; isMA; false]);
start_idx = find(d == 1);         % MA starts
end_idx   = find(d == -1) - 1;    % MA ends

nMA = numel(start_idx);
ma_on_s  = nan(nMA,1);
ma_off_s = nan(nMA,1);
state_before = strings(nMA,1);

for k = 1:nMA
    i0 = start_idx(k);
    i1 = end_idx(k);

    % onset/offset time in seconds (0-based)
    ma_on_s(k)  = (i0-1) * epoch_sec;
    ma_off_s(k) = (i1)   * epoch_sec;   % end = last sample * epoch_sec

    % ---- state_before: last non-MA sample before i0 ----
    pre_idx = i0 - 1;
    while pre_idx >= 1 && codes_vec(pre_idx) == codes.MA
        pre_idx = pre_idx - 1;
    end

    if pre_idx < 1
        state_before(k) = "UNK";
    else
        c = codes_vec(pre_idx);
        if isfield(codes,'NREM') && c == codes.NREM
            state_before(k) = "NREM";
        elseif isfield(codes,'REM') && c == codes.REM
            state_before(k) = "REM";
        elseif isfield(codes,'WK') && c == codes.WK
            state_before(k) = "WK";
        else
            state_before(k) = "OTHER";
        end
    end
end

% --------- Wrap into table ----------
mouse_id  = string(mouse_id);
genotype  = string(genotype);
condition = string(condition);

MA = table( ...
    repmat(mouse_id,  nMA,1), ...
    repmat(genotype,  nMA,1), ...
    repmat(condition, nMA,1), ...
    ma_on_s, ...
    ma_off_s, ...
    state_before, ...
    'VariableNames',{'mouse','genotype','condition','ma_on_s','ma_off_s','state_before'});
end
