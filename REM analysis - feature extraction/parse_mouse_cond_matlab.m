function [mouse, cond] = parse_mouse_cond_matlab(fname)
s = lower(string(fname));
mouse = "unknownMouse"; cond = "unknownCond";
m = regexp(s,'(mouse\d+)','tokens','once'); if ~isempty(m), mouse=string(m{1}); end
tokens = regexp(s,'[\w]+','match'); toks = string(unique(tokens));
if any(ismember(toks,["baseline","base"])), cond="baseline"; return; end
if any(ismember(toks,["ambtemp","ambient","amb_temp","roomtemp","rt"])), cond="ambtemp"; return; end
if any(ismember(toks,["drugs","drug","cno","ket","muscimol"])), cond="drugs"; return; end
% fallback: token after mouseXX separated by - or _
if ~strcmp(mouse,"unknownMouse")
    tail = extractAfter(s, mouse);
    m2 = regexp(tail,'[-_]*([a-z0-9]+)','tokens','once');
    if ~isempty(m2), cond = string(m2{1}); end
end
end
