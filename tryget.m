function val = tryget(S, keys)
val = [];
for i=1:numel(keys)
    if isfield(S, keys{i})
        v = S.(keys{i});
        if isnumeric(v) && isscalar(v)
            val = v; return;
        end
    end
end
end

