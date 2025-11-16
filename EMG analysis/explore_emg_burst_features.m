function explore_emg_burst_features(EMG_GROUP)
% EXPLORE_EMG_BURST_FEATURES
% -------------------------------------------------------------------------
% Visualise EMG burst timing and amplitude in baseline vs ambtemp
% to help separate physiology from artifact.
%
% INPUT
%   EMG_GROUP : struct from run_emg_batch_auto

sessions = EMG_GROUP.sessions;
out      = EMG_GROUP.out;

% collect bursts for baseline vs ambtemp
IBI_baseline = [];
AMP_baseline = [];

IBI_ambtemp  = [];
AMP_ambtemp  = [];   % aligned so IBI_ambtemp(k) ↔ AMP_ambtemp(k)

for k = 1:numel(sessions)
    s   = sessions(k);
    OUT = out{k};
    if ~isfield(OUT,'bursts') || isempty(OUT.bursts)
        continue;
    end
    b = OUT.bursts;

    % mid times and IBI
    mid = (b.start_s + b.end_s)/2;
    if numel(mid) < 2, continue; end
    ibi = diff(mid);
    amp = b.peak_amp;

    if strcmpi(s.cond,'baseline')
        % for baseline, we just want overall distributions, so no need to align
        IBI_baseline = [IBI_baseline; ibi(:)];   %#ok<AGROW>
        AMP_baseline = [AMP_baseline; amp(:)];   %#ok<AGROW>
    elseif strcmpi(s.cond,'ambtemp')
        % for ambtemp, we want aligned IBI↔amp pairs:
        % IBI(i) is between burst i and i+1, so pair with amplitude of burst i+1
        IBI_ambtemp = [IBI_ambtemp; ibi(:)];          %#ok<AGROW>
        AMP_ambtemp = [AMP_ambtemp; amp(2:end)];      %#ok<AGROW>
    end
end

figure('Name','EMG burst features: baseline vs ambtemp','Color','w');

% IBI histograms
subplot(2,2,1);
hold on;
if ~isempty(IBI_baseline)
    histogram(IBI_baseline, 0:0.05:5, 'Normalization','pdf', ...
        'FaceAlpha',0.5,'DisplayName','baseline');
end
if ~isempty(IBI_ambtemp)
    histogram(IBI_ambtemp, 0:0.05:5, 'Normalization','pdf', ...
        'FaceAlpha',0.5,'DisplayName','ambtemp');
end
xlabel('IBI (s)');
ylabel('PDF');
title('Inter-burst interval distribution (0–5 s)');
legend; box off;

% zoom on short IBIs
subplot(2,2,2);
hold on;
if ~isempty(IBI_baseline)
    histogram(IBI_baseline, 0:0.02:1, 'Normalization','pdf', ...
        'FaceAlpha',0.5,'DisplayName','baseline');
end
if ~isempty(IBI_ambtemp)
    histogram(IBI_ambtemp, 0:0.02:1, 'Normalization','pdf', ...
        'FaceAlpha',0.5,'DisplayName','ambtemp');
end
xlabel('IBI (s)');
ylabel('PDF');
title('Inter-burst interval (0–1 s)');
legend; box off;

% amplitude distributions
subplot(2,2,3);
hold on;
if ~isempty(AMP_baseline)
    histogram(AMP_baseline, 50, 'Normalization','pdf', ...
        'FaceAlpha',0.5,'DisplayName','baseline');
end
if ~isempty(AMP_ambtemp)
    histogram(AMP_ambtemp, 50, 'Normalization','pdf', ...
        'FaceAlpha',0.5,'DisplayName','ambtemp');
end
xlabel('Burst peak amplitude (envelope units)');
ylabel('PDF');
title('Burst amplitude');
legend; box off;

% IBI vs amplitude scatter (ambtemp)
subplot(2,2,4);
if ~isempty(IBI_ambtemp) && ~isempty(AMP_ambtemp)
    scatter(IBI_ambtemp, AMP_ambtemp, 5, 'k','filled');
    xlabel('IBI (s)');
    ylabel('Peak amplitude');
    title('Ambtemp bursts: IBI vs amplitude');
    box off;
else
    text(0.5,0.5,'No ambtemp bursts found','HorizontalAlignment','center');
    axis off;
end

sgtitle('EMG burst features in baseline vs ambtemp');

end
