EMG_GROUP = run_emg_batch_auto( ...
    '/Users/margaridaseabra/signalnotscored', ...
    '/Users/margaridaseabra/15.11scores');
%%
run_emg_group_plots(EMG_GROUP);

%%
explore_emg_burst_features(EMG_GROUP);
%%
idxAmb = find(strcmpi({EMG_GROUP.sessions.cond}, 'ambtemp'));
k = idxAmb(1);   % choose whichever you want to inspect
%%
OUT = EMG_GROUP.out{k};  % some ambtemp session
b   = OUT.bursts;

mid = (b.start_s + b.end_s)/2;
IBI = diff(mid);

% rough IBI mode estimation for this session:
edges = 0:0.05:5;
[counts,edges] = histcounts(IBI, edges);
[~,mx] = max(counts);
ibi_mode = mean(edges(mx:mx+1));  % center of most common bin

artifact_idx = flag_emg_artifact_bursts(b, ibi_mode, 0.05, 0.2);

% now you can see how many:
fprintf('Session: %d bursts, %d flagged as artifact (%.1f%%)\n', ...
        height(b), nnz(artifact_idx), 100*nnz(artifact_idx)/height(b));

% if you like what you see, you can drop them:
b_clean = b(~artifact_idx,:);
