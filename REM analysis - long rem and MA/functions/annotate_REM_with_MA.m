function REM = annotate_REM_with_MA(BOUTS, MA_tbl, varargin)
% annotate_REM_with_MA
% -------------------------------------------------------------
% Keep only REM bouts and annotate:
%   - is_long (duration >= quantile within recording)
%   - n_MA_pre, MA_dur_pre within pre-REM window
%
% Inputs:
%   BOUTS  : table from build_bout_table
%   MA_tbl : table with at least start_s, end_s (in seconds)
%
% Optional name-value:
%   'pre_window_s'  : e.g. 60
%   'long_quantile' : e.g. 0.75

ip = inputParser;
addParameter(ip,'pre_window_s',60,@isscalar);
addParameter(ip,'long_quantile',0.75,@(x)x>0 && x<1);
parse(ip,varargin{:});

pre_window_s  = ip.Results.pre_window_s;
long_q        = ip.Results.long_quantile;

% keep REM only
isREM = (BOUTS.state=="REM");
REM = BOUTS(isREM,:);
if isempty(REM); return; end

% --- define "long REM" within this recording ----------------------------
thr = quantile(REM.dur_s, long_q);
REM.is_long = REM.dur_s >= thr;

% --- pre-REM MA counts ---------------------------------------------------
if ~isempty(MA_tbl)
    if ~all(ismember({'start_s','end_s'}, MA_tbl.Properties.VariableNames))
        error('MA_tbl must have columns start_s and end_s (seconds).');
    end
    MA_start = MA_tbl.start_s(:);
    MA_end   = MA_tbl.end_s(:);
else
    MA_start = [];
    MA_end   = [];
end

nB = height(REM);
n_MA_pre   = zeros(nB,1);
MA_dur_pre = zeros(nB,1);

for iB = 1:nB
    tREM = REM.t_start(iB);
    win_start = tREM - pre_window_s;
    win_end   = tREM;
    
    if isempty(MA_start)
        n_MA_pre(iB)   = 0;
        MA_dur_pre(iB) = 0;
        continue;
    end
    
    inwin = (MA_start >= win_start) & (MA_start < win_end);
    n_MA_pre(iB)   = sum(inwin);
    
    if any(inwin)
        MA_dur_pre(iB) = sum(MA_end(inwin) - MA_start(inwin));
    else
        MA_dur_pre(iB) = 0;
    end
end

REM.n_MA_pre   = n_MA_pre;
REM.MA_dur_pre = MA_dur_pre;
REM.long_thr_s = repmat(thr, nB, 1);   % keep threshold for reference
