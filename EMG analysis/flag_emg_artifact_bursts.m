function artifact_idx = flag_emg_artifact_bursts(bursts, ibi_mode, ibi_tol, amp_cv_max)
% FLAG_EMG_ARTIFACT_BURSTS
% -------------------------------------------------------------------------
% Simple rule-based artifact detector for EMG bursts:
%   - marks bursts as artifact if their IBI is near a given mode
%     (ibi_mode ± ibi_tol), and their amplitude distribution is "tight".
%
% INPUTS
%   bursts     : table OUT.bursts from emg_burst_analysis
%   ibi_mode   : dominant inter-burst interval (s) to treat as artifact
%   ibi_tol    : tolerance around ibi_mode (e.g. 0.05 s)
%   amp_cv_max : max coefficient of variation of amplitudes to consider the
%                cluster "stereotyped" enough (e.g. 0.2)
%
% OUTPUT
%   artifact_idx : logical vector (height(bursts) x 1), true = artifact

nB = height(bursts);
artifact_idx = false(nB,1);
if nB < 3
    return;
end

mid = (bursts.start_s + bursts.end_s)/2;
IBI = diff(mid);

% bursts have IBI defined between them: associate IBI with the *later* burst
ibi_flag = false(nB,1);
ibi_flag(2:end) = abs(IBI - ibi_mode) < ibi_tol;

% check amplitude stability within this group
amp = bursts.peak_amp;
amp_art = amp(ibi_flag);
if numel(amp_art) < 5
    % not enough candidate artifacts, skip
    return;
end

cv = std(amp_art,'omitnan') / max(mean(amp_art,'omitnan'), eps);
if cv > amp_cv_max
    % amplitudes too variable → probably not pure artifact
    return;
end

artifact_idx = ibi_flag;
end
