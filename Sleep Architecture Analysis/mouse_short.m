function lbl = mouse_short(s)
% MOUSE_SHORT  Make compact labels like "m1" from "mouse1"/"Mouse-03"/"03".
%
% INPUT
%   s : string / char / cellstr with mouse IDs
%
% OUTPUT
%   lbl : same size as s, labels like "m1", "m03", etc.

s = string(s);
lbl = strings(size(s));

for k = 1:numel(s)
    t = s(k);

    % Get trailing digits (e.g. "mouse10" -> "10")
    num = regexp(t, '\d+$', 'match', 'once');

    if ~isempty(num)
        % Use "m" + those digits
        lbl(k) = "m" + string(num);
    else
        % Fallback: trimmed original, max 8 chars
        t = regexprep(t, '^\s+|\s+$', '');  % trim
        if strlength(t) > 8
            t = extractBefore(t, 9);
        end
        lbl(k) = t;
    end
end
end
