function [st,en] = runs_from_codes(code)
% Return start/end indices of contiguous runs in integer vector "code"
code = code(:);
if isempty(code)
    st = []; en = [];
    return;
end
d = diff(code);
change_idx = [1; find(d ~= 0) + 1];
st = change_idx;
en = [change_idx(2:end)-1; numel(code)];
end
