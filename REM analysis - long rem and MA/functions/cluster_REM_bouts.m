function CLUST = cluster_REM_bouts(REM, varargin)
% cluster_REM_bouts
% -------------------------------------------------------------
% Cluster REM bouts based on inter-REM interval and classify clusters as:
%   - ShortOnly  : only normal REMs, no long ones
%   - ShortLong  : >=1 short REM followed by long REM
%   - LongOnly   : cluster with single long REM only
%
% Inputs:
%   REM : table from annotate_REM_with_MA (must include t_start, t_end,
%         dur_s, is_long)
%
% Optional:
%   'max_gap_s' : max gap between bouts to be in same cluster (default 60)

ip = inputParser;
addParameter(ip,'max_gap_s',60,@isscalar);
parse(ip,varargin{:});
max_gap_s = ip.Results.max_gap_s;

CLUST = table();
if isempty(REM); return; end

% sort by start time (just in case)
REM = sortrows(REM, 't_start');

nB = height(REM);
cluster_id = zeros(nB,1);
cur_id = 1;
cluster_id(1) = cur_id;

for iB = 2:nB
    gap = REM.t_start(iB) - REM.t_end(iB-1);
    if gap <= max_gap_s
        cluster_id(iB) = cur_id;
    else
        cur_id = cur_id + 1;
        cluster_id(iB) = cur_id;
    end
end

REM.cluster_id = cluster_id;

% now summarise per cluster
unique_ids = unique(cluster_id);
nC = numel(unique_ids);

mouse_col  = strings(nC,1);
geno_col   = strings(nC,1);
cond_col   = strings(nC,1);
cluster_type = strings(nC,1);
n_REM_cluster = zeros(nC,1);
n_short_before_long = nan(nC,1);
cluster_t_start = nan(nC,1);
cluster_t_end   = nan(nC,1);
cluster_total_dur = nan(nC,1);

for iC = 1:nC
    cid = unique_ids(iC);
    idx = (cluster_id == cid);
    R = REM(idx,:);
    
    mouse_col(iC) = R.mouse(1);
    geno_col(iC)  = R.geno(1);
    cond_col(iC)  = R.cond(1);
    
    n_REM_cluster(iC) = height(R);
    cluster_t_start(iC) = min(R.t_start);
    cluster_t_end(iC)   = max(R.t_end);
    cluster_total_dur(iC) = sum(R.dur_s);
    
    has_long = any(R.is_long);
    
    if has_long
        if height(R)==1  % only one long REM
            cluster_type(iC) = "LongOnly";
            n_short_before_long(iC) = 0;
        else
            % find first long REM and count short before it
            first_long_idx = find(R.is_long, 1, 'first');
            n_short_before_long(iC) = sum(~R.is_long(1:first_long_idx-1));
            cluster_type(iC) = "ShortLong";
        end
    else
        cluster_type(iC) = "ShortOnly";
        n_short_before_long(iC) = NaN;
    end
end

CLUST = table;
CLUST.mouse  = mouse_col;
CLUST.geno   = geno_col;
CLUST.cond   = cond_col;
CLUST.cluster_id = unique_ids;
CLUST.cluster_type = cluster_type;
CLUST.n_REM_cluster = n_REM_cluster;
CLUST.n_short_before_long = n_short_before_long;
CLUST.t_start = cluster_t_start;
CLUST.t_end   = cluster_t_end;
CLUST.total_REM_dur_s = cluster_total_dur;
