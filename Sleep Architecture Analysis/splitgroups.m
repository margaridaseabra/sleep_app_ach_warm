function [cellOut, varargout] = splitgroups(varargin)
% Utility: reconstruct group keys from grouping vector
G = varargin{end};
K = numel(varargin)-1;
keys = cell(1,K);
for k=1:K
    v = varargin{k};
    if iscell(v), v = string(v); end
    keys{k} = splitapply(@(x){x(1)}, v, G);
end
cellOut = keys{1};
for k=2:K, varargout{k-1} = keys{k}; end
end
