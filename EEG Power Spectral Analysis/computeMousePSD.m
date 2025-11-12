function result = computeMousePSD(eegFile, scoreFile, genotype, mouseID)
% Computes averaged PSD per sleep state for a single mouse

% Load EEG and scoring data
eegData = load(eegFile);
scoreData = load(scoreFile);

EEG = eegData.eeg;
Fs  = eegData.eeg_frequency;

TT = scoreData.TT;
C  = scoreData.C;

score = TT{:,1};
score(score == C.MA) = C.NREM; % Merge MA into NREM

epochLength = scoreData.epochSec;
samplesPerEpoch = round(Fs * epochLength);

states = struct('Wake', C.WK, 'NREM', C.NREM, 'REM', C.REM);
stateNames = fieldnames(states);

result.mouseID = mouseID;
result.genotype = genotype;

for s = 1:length(stateNames)
    stateName = stateNames{s};
    stateVal = states.(stateName);
    epochs = find(score == stateVal);
    
    allPxx = [];
    F_ref = [];

    for i = 1:length(epochs)
        idx1 = round((epochs(i)-1)*samplesPerEpoch + 1);
        idx2 = idx1 + samplesPerEpoch - 1;
        if idx2 <= length(EEG)
            segment = EEG(idx1:idx2);
            [pxx, F] = pwelch(segment, hamming(512), [], [], Fs);
            if isempty(F_ref), F_ref = F; end
            if length(pxx) == length(F_ref)
                allPxx(:, end+1) = pxx;
            end
        end
    end

    result.(stateName).F = F_ref;
    result.(stateName).Pxx = mean(allPxx, 2);
end
end
