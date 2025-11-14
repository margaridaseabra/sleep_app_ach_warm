function s = safe_str(s)
if isempty(s), s = 'NA'; return; end
s = regexprep(s,'\s+','_');
s = regexprep(s,'[^\w-]','');
end
