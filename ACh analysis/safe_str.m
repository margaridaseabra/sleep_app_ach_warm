function s = safe_str(s)
% make a string safe for filenames
if isempty(s), s = 'NA'; return; end
s = regexprep(s,'\s+','_');
s = regexprep(s,'[^\w-]','');
end

