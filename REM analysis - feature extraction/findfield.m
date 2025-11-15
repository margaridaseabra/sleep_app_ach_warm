function name = findfield(S, candidates)
for i=1:numel(candidates)
    if isfield(S, candidates{i}), name = candidates{i}; return; end
end
error('None of the candidate fields found: %s', strjoin(candidates,', '));
end

function name = tryfield(S, candidates)
name = '';
for i=1:numel(candidates)
    if isfield(S, candidates{i}), name = candidates{i}; return; end
end
end

