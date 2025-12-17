function [REM, edges_s] = add_REM_size_category(REM, varargin)
% add_REM_size_category  Add Short / Medium / Long REM size category.
%
% Usage:
%   [REM, edges_s] = add_REM_size_category(REM)
%   [REM, edges_s] = add_REM_size_category(REM,'durVar','dur_s')
%   [REM, edges_s] = add_REM_size_category(REM,'edges',[20 60])
%
% Inputs:
%   REM   : table with a column containing REM bout duration (seconds)
%
% Options:
%   'durVar' : name of the duration column in REM (char or string).
%              If omitted, tries dur_s, dur, REM_dur_s in this order.
%   'edges'  : [short_medium_boundary, medium_long_boundary] in seconds.
%              If omitted, tertiles of the pooled REM durations are used.
%
% Outputs:
%   REM      : same table with new categorical variable REM_size
%              ('Short','Medium','Long')
%   edges_s  : 1x2 vector with the thresholds actually used (s)

p = inputParser;
addParameter(p,'durVar','',@(x)ischar(x)||isstring(x));
addParameter(p,'edges',[],@(x)isnumeric(x)&&numel(x)==2);
parse(p,varargin{:});
durVar  = string(p.Results.durVar);
edges_s = p.Results.edges;

% ---------- figure out which column holds duration ----------
if durVar ~= ""
    % User-specified column
    if ~ismember(durVar, string(REM.Properties.VariableNames))
        error('add_REM_size_category:BadDurVar', ...
            'Table does not have a column named "%s". Available vars: %s', ...
            durVar, strjoin(REM.Properties.VariableNames, ', '));
    end
else
    % Try common names automatically
    candidates = ["dur_s","dur","REM_dur_s","bout_dur_s"];
    found = candidates(ismember(candidates, string(REM.Properties.VariableNames)));
    if isempty(found)
        error('add_REM_size_category:NoDurColumn', ...
            ['Could not find a duration column. Tried: %s\nAvailable vars: %s'], ...
            strjoin(candidates, ', '), ...
            strjoin(REM.Properties.VariableNames, ', '));
    end
    durVar = found(1);  % first match
end

dur = REM.(durVar);

% ---------- compute thresholds ----------
if isempty(edges_s)
    edges_s = quantile(dur, [1/3 2/3]);
end

% ---------- assign Short / Medium / Long ----------
sz = repmat("Medium", height(REM), 1);
sz(dur <  edges_s(1)) = "Short";
sz(dur >  edges_s(2)) = "Long";

REM.REM_size = categorical(sz, ["Short","Medium","Long"]);

end
