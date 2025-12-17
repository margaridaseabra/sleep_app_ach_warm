function H = ach_group_state_heatmap(groupTransMat, stateName, geno, cond, t_pre, t_post, dt, zscore_rows)
% ACH_GROUP_STATE_HEATMAP
%   Group heatmap for a given state, genotype and condition.
%
%   H = ach_group_state_heatmap('ACh_periTransitions_all.mat',...
%                               'NREM','WT','drugs');

if nargin < 1 || isempty(groupTransMat)
    [fname,fpath] = uigetfile('*.mat','Select ACh_periTransitions_all.mat');
    if isequal(fname,0), error('No file selected.'); end
    groupTransMat = fullfile(fpath,fname);
end
if nargin < 2 || isempty(stateName), stateName = 'NREM';    end
if nargin < 3 || isempty(geno),       geno      = 'WT';     end
if nargin < 4 || isempty(cond),       cond      = 'drugs';  end
if nargin < 5 || isempty(t_pre),      t_pre     = 30;       end
if nargin < 6 || isempty(t_post),     t_post    = 60;       end
if nargin < 7 || isempty(dt),         dt        = 0.1;      end
if nargin < 8 || isempty(zscore_rows),zscore_rows = true;   end

S = load(groupTransMat,'GROUP_TRANS');
GROUP_TRANS = S.GROUP_TRANS;

sessions = GROUP_TRANS.sessions;
OUT_all  = GROUP_TRANS.out;
nSess    = numel(sessions);

t_common = -t_pre:dt:t_post;
allEvents = [];
usedSess  = false(nSess,1);

for k = 1:nSess
    s = sessions(k);
    if ~strcmp(s.geno, geno) || ~strcmp(s.cond, cond)
        continue;
    end

    OUT_T = OUT_all{k};
    if isempty(OUT_T) || ~isfield(OUT_T,stateName), continue; end
    M = OUT_T.(stateName).M;
    if isempty(M), continue; end

    t_sess = OUT_T.t_rel(:).';
    if numel(t_sess) ~= size(M,2), continue; end

    % interpolate each event onto common time axis
    M_interp = nan(size(M,1), numel(t_common));
    for e = 1:size(M,1)
        M_interp(e,:) = interp1(t_sess, M(e,:), t_common, 'linear', NaN);
    end

    allEvents = [allEvents; M_interp]; %#ok<AGROW>
    usedSess(k) = true;
end

if isempty(allEvents)
    warning('No events for %s onset in %s, %s.', stateName, geno, cond);
    H = [];
    return;
end

if zscore_rows
    mu = mean(allEvents,2,'omitnan');
    sd = std(allEvents,[],2,'omitnan');
    allEvents = (allEvents - mu) ./ sd;
end

H.M     = allEvents;
H.t     = t_common;
H.Nev   = size(allEvents,1);
H.nSess = sum(usedSess);

% ---- plot ----
figure('Color','w','Position',[100 100 650 500]);
imagesc(t_common, 1:H.Nev, allEvents);
axis xy;
xlabel('Time from onset (s)');
ylabel('Event #');
title(sprintf('%s onset, %s %s (events from %d sessions, n=%d)', ...
    stateName, geno, cond, H.nSess, H.Nev));
colormap hot; colorbar;
hold on;
plot([0 0],[0 H.Nev+1],'w--','LineWidth',1.2,'HandleVisibility','off');

end
