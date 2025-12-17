function bouts = build_context_bouts(score, t_scores, codes)
% Build bouts with context info (next and previous states).
%
% OUTPUT: bouts(k) with fields:
%   .state        : 'Wake'/'NREM'/'REM'/'MA'
%   .code         : numeric code
%   .start_s, .end_s
%   .prev_state   : 'Wake'/'NREM'/'REM'/'MA'/'none'
%   .next_state   : same
%
% Requires t_scores in seconds, same length as score.

stateMap = containers.Map( ...
    [codes.WK, codes.NREM, codes.REM, codes.MA], ...
    {'Wake','NREM','REM','MA'});

dt = median(diff(t_scores));
score = score(:);
t_scores = t_scores(:);

% find contiguous runs
d  = diff([0; score; 0] ~= 0);  % 1 on change, -1 on end
i1 = find(d ==  1);
i2 = find(d == -1)-1;

nb = numel(i1);
bouts = struct('state',[],'code',[],'start_s',[],'end_s',[], ...
               'prev_state',[],'next_state',[]);

for k = 1:nb
    code = score(i1(k));
    if ~isKey(stateMap, code), continue; end
    bouts(k).code    = code;
    bouts(k).state   = stateMap(code);
    bouts(k).start_s = t_scores(i1(k));
    bouts(k).end_s   = t_scores(i2(k)) + dt;
end

% context tags (skip MA when looking for neighbours)
for k = 1:nb
    % look backwards for previous non-MA
    prev = 'none';
    for j = k-1:-1:1
        if bouts(j).code ~= codes.MA
            prev = bouts(j).state;
            break;
        end
    end
    bouts(k).prev_state = prev;

    % look forwards for next non-MA
    nxt = 'none';
    for j = k+1:nb
        if bouts(j).code ~= codes.MA
            nxt = bouts(j).state;
            break;
        end
    end
    bouts(k).next_state = nxt;
end
end
