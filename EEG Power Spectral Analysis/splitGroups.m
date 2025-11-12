function groups = splitGroups(indices)
    if isempty(indices)
        groups = {};
        return;
    end
    d = diff(indices);
    boundaries = [0; find(d > 1); length(indices)];
    groups = cell(length(boundaries)-1, 1);
    for i = 1:length(groups)
        groups{i} = indices((boundaries(i)+1):boundaries(i+1));
    end
end