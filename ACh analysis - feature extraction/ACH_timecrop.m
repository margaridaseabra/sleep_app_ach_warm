function seg = ACH_timecrop(ACH, fs_ach, fs_eeg, idx_eeg)
if isempty(ACH) || isempty(fs_ach) || isempty(fs_eeg)
    seg = [];
    return;
end
t = (idx_eeg-1)/fs_eeg;
idx_ach = round(t*fs_ach) + 1;
idx_ach(idx_ach<1) = 1;
idx_ach(idx_ach>numel(ACH)) = numel(ACH);
seg = ACH(idx_ach);
end