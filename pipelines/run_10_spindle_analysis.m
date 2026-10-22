%% RUN_10_SPINDLE_ANALYSIS
% Publication-oriented spindle analysis for all baseline, warming and
% donepezil sessions in the APP/PS1 ambient-temperature project.
%
% Reference detector:
%   Bandarabadi et al., Nature Communications (2020)
%   "A role for spindles in the onset of rapid eye movement sleep"
%
% PRIMARY detector choices
% ------------------------
%   EEG analysis Fs                200 Hz
%   mains notch                    50 Hz
%   wavelet                        fbsp2-1-2
%   wavelet-energy band            9-16 Hz, 0.5-Hz centers
%   extra 1/f weighting            NO (paper-primary analysis)
%   energy smoothing               200-ms Hann
%   low threshold                  mean + 1 SD
%   high threshold                 mean + 3 SD
%   rejection threshold            mean + 20 SD
%   threshold outlier exclusion    > mean + 10 SD
%   duration                       0.4-2.0 s
%   NREM-like overlap              >= 80%
%   cycles                         5-30
%   frequency specificity          spindle-band peak must exceed peaks in
%                                  6-8.5 and 16.5-20 Hz
%
% NREM-like definition:
%   score 1 (NREM) + score 15 (MA).
%   Before analysis, Wake bouts <=15 s are reclassified to MA, matching
%   the sleep-architecture pipeline.
%
% PRIMARY transition endpoint:
%   spindle rate during the final 25 s before NREM->REM vs NREM->Wake.
%   A descriptive 60-s, 1-Hz occupancy trajectory is also exported.
%
% Output:
%   outputs/eeg/spindles/
%       spindle_events.csv
%       spindle_window_summary.csv
%       spindle_transition_bouts.csv
%       spindle_transition_mouse_summary.csv
%       spindle_transition_timecourse_bouts.csv
%       spindle_transition_timecourse_mouse.csv
%       spindle_session_qc.csv
%       qc_figures/
%
% Run:
%   cd('/path/to/sleep_app_ach_warm')
%   run('pipelines/run_10_spindle_analysis.m')

clear;
clc;

fprintf('\n============================================================\n');
fprintf('SPINDLE ANALYSIS\n');
fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
% Setup
% -------------------------------------------------------------------------

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);

addpath(fullfile(repoDir,'config'));
P = paths_config();

sessionsFile = P.sessions;
miceFile = P.mice;

standardDir = fullfile(P.sleep_output,'all','standardized_scores');

assert(isfile(sessionsFile),'Missing %s',sessionsFile);
assert(isfile(miceFile),'Missing %s',miceFile);
assert(isfolder(standardDir), ...
    'Missing standardized scores. Run run_01_sleep_architecture.m first.');

outDir = fullfile(P.eeg_output,'spindles');
figDir = fullfile(outDir,'qc_figures');

if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(figDir), mkdir(figDir); end

warmingTimesFile = fullfile(repoDir,'config','warming_times.csv');
drugTimesFile = fullfile(repoDir,'config','drug_times.csv');

%% ------------------------------------------------------------------------
% Locked parameters
% -------------------------------------------------------------------------

CFG.target_fs = 200;

CFG.notch_hz = 50;
CFG.notch_Q = 35;

CFG.wavelet = 'fbsp2-1-2';
CFG.wavelet_freqs_hz = 9:0.5:16;
CFG.weight_wavelet_by_frequency = false;

CFG.energy_smooth_s = 0.200;

CFG.threshold_low_sd = 1;
CFG.threshold_high_sd = 3;
CFG.threshold_max_sd = 20;
CFG.threshold_outlier_sd = 10;

CFG.min_duration_s = 0.4;
CFG.max_duration_s = 2.0;
CFG.min_nrem_like_fraction = 0.80;

CFG.min_cycles = 5;
CFG.max_cycles = 30;

CFG.spindle_band_hz = [9 16];
CFG.low_flank_hz = [6 8.5];
CFG.high_flank_hz = [16.5 20];

% Chunk the CWT so a 9-h recording never creates one giant coefficient
% matrix in memory. Margin reduces edge effects at chunk boundaries.
CFG.cwt_chunk_s = 300;
CFG.cwt_margin_s = 5;

CFG.transition_primary_s = 25;
CFG.transition_timecourse_s = 60;

CFG.ma_threshold_s = 15;

CFG.save_qc_plots = true;
CFG.max_qc_plots = 12;

%% ------------------------------------------------------------------------
% Metadata
% -------------------------------------------------------------------------

Sessions = readtable(sessionsFile, ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

Mice = readtable(miceFile, ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

Sessions = Sessions(to_logical(Sessions.include_sleep),:);

WARM = table();
DRUG = table();

if isfile(warmingTimesFile)
    WARM = readtable(warmingTimesFile, ...
        'TextType','string', ...
        'VariableNamingRule','preserve');
end

if isfile(drugTimesFile)
    DRUG = readtable(drugTimesFile, ...
        'TextType','string', ...
        'VariableNamingRule','preserve');
end

%% ------------------------------------------------------------------------
% Collectors
% -------------------------------------------------------------------------

EVENTS = table();
WINDOWS = table();
TRANSITIONS = table();
TRANSITION_TC = table();
QC = table();

qcPlotCount = 0;

%% ------------------------------------------------------------------------
% Session loop
% -------------------------------------------------------------------------

for i = 1:height(Sessions)

    sid = string(Sessions.session_id(i));
    mouseID = string(Sessions.mouse_id(i));
    condition = string(Sessions.condition(i));
    cohort = string(Sessions.cohort(i));

    mrow = find(string(Mice.mouse_id)==mouseID,1);

    if isempty(mrow)
        error('No mouse metadata for %s',mouseID);
    end

    genotype = string(Mice.genotype(mrow));

    sex = missing;
    if ismember("sex",string(Mice.Properties.VariableNames))
        sex = string(Mice.sex(mrow));
    end

    cohortDisplay = cohort;
    if cohortDisplay=="thesis"
        cohortDisplay="original";
    end

    fprintf('[%02d/%02d] %-31s %-10s %-4s ', ...
        i,height(Sessions),sid,condition,genotype);

    %% Resolve raw EEG recording

    rawRel = "";

    if ismember("recording_file",string(Sessions.Properties.VariableNames))
        rawRel = string(Sessions.recording_file(i));
    end

    if ismissing(rawRel) || strlength(strtrim(rawRel))==0 || ...
            ~isfile(fullfile(P.data_root,char(rawRel)))

        rawRel = infer_raw_recording_path( ...
            Sessions(i,:),P.data_root);

    end

    rawFile = fullfile(P.data_root,char(rawRel));

    assert(isfile(rawFile), ...
        'Raw EEG recording not found for %s:\n%s',sid,rawFile);

    scoreFile = fullfile(standardDir,char(sid + ".csv"));

    assert(isfile(scoreFile), ...
        'Standardized score file not found:\n%s',scoreFile);

    %% Load EEG

    vars = string({whos('-file',rawFile).name});

    assert(ismember("eeg",vars), ...
        'No top-level eeg variable in %s',rawFile);

    assert(ismember("eeg_frequency",vars), ...
        'No top-level eeg_frequency variable in %s',rawFile);

    R = load(rawFile,'eeg','eeg_frequency');

    eeg = double(R.eeg(:));
    rawFs = double(R.eeg_frequency);

    assert(isscalar(rawFs) && isfinite(rawFs) && rawFs>0, ...
        'Invalid eeg_frequency in %s',rawFile);

    %% Load 1-Hz scoring

    ST = readtable(scoreFile, ...
        'TextType','string', ...
        'VariableNamingRule','preserve');

    assert(all(ismember(["time_s","score"], ...
        string(ST.Properties.VariableNames))), ...
        'Standardized scores require time_s and score.');

    scoreTime = double(ST.time_s(:));
    scores = double(ST.score(:));

    keep = isfinite(scoreTime) & isfinite(scores);

    scoreTime = scoreTime(keep);
    scores = scores(keep);

    assert(~isempty(scores),'No valid sleep scores for %s',sid);

    scoreTime = scoreTime-scoreTime(1);

    % Match sleep pipeline: short Wake <=15 s becomes MA.
    scores = reclassify_short_wake(scores,CFG.ma_threshold_s);

    %% EEG preprocessing

    if abs(rawFs-CFG.target_fs)>1e-9
        [pRat,qRat] = rat(CFG.target_fs/rawFs,1e-10);
        eeg = resample(eeg,pRat,qRat);
    end

    fs = CFG.target_fs;
    eeg = double(eeg(:));
    eeg = detrend(eeg,'constant');

    if CFG.notch_hz < fs/2
        eeg = apply_notch_compat( ...
            eeg,fs,CFG.notch_hz,CFG.notch_Q);
    end

    N = numel(eeg);
    durationS = N/fs;

    %% Map scores to EEG samples

    [sampleState,sampleValid] = scores_to_samples( ...
        scoreTime,scores,N,fs);

    nremLike = sampleValid & ...
        (sampleState==1 | sampleState==15);

    nremLikeS = sum(nremLike)/fs;

    if nremLikeS < 30
        warning('Only %.1f s NREM-like data in %s; skipping.', ...
            nremLikeS,sid);
        continue;
    end

    %% Wavelet energy

    energy = fbsp_energy_chunked( ...
        eeg,fs,CFG.wavelet,CFG.wavelet_freqs_hz, ...
        CFG.weight_wavelet_by_frequency, ...
        CFG.cwt_chunk_s,CFG.cwt_margin_s);

    nSmooth = max(3,round(CFG.energy_smooth_s*fs));

    if mod(nSmooth,2)==0
        nSmooth=nSmooth+1;
    end

    h = hann(nSmooth);
    h = h/sum(h);

    energy = conv(energy,h,'same');

    %% Thresholds from NREM-like samples

    eThr = energy(nremLike & isfinite(energy));

    mu0 = mean(eThr);
    sd0 = std(eThr);

    if ~isfinite(sd0) || sd0<=0
        warning('Invalid spindle-energy SD for %s; skipping.',sid);
        continue;
    end

    eClean = eThr(eThr <= mu0 + CFG.threshold_outlier_sd*sd0);

    if numel(eClean)<fs*30
        eClean=eThr;
    end

    muE = mean(eClean);
    sdE = std(eClean);

    thrL = muE + CFG.threshold_low_sd*sdE;
    thrH = muE + CFG.threshold_high_sd*sdE;
    thrM = muE + CFG.threshold_max_sd*sdE;

    %% Bandpassed spindle signal

    filtOrder = max(20,round(3*fs/CFG.spindle_band_hz(1)));

    if mod(filtOrder,2)==1
        filtOrder=filtOrder+1;
    end

    bSp = fir1( ...
        filtOrder, ...
        CFG.spindle_band_hz/(fs/2), ...
        'bandpass', ...
        hann(filtOrder+1));

    eegSp = filtfilt(bSp,1,eeg);

    %% Detect events

    candidateRuns = logical_runs(energy>thrH);

    expanded = zeros(0,2);

    for r = 1:size(candidateRuns,1)

        a = candidateRuns(r,1);
        b = candidateRuns(r,2);

        while a>1 && energy(a-1)>thrL
            a=a-1;
        end

        while b<N && energy(b+1)>thrL
            b=b+1;
        end

        expanded(end+1,:)=[a b]; %#ok<AGROW>
    end

    expanded = merge_intervals(expanded);

    sessionEvents = table();

    for r = 1:size(expanded,1)

        a = expanded(r,1);
        b = expanded(r,2);

        dur = (b-a+1)/fs;

        if dur<CFG.min_duration_s || dur>CFG.max_duration_s
            continue;
        end

        eSeg = energy(a:b);

        if max(eSeg)>thrM
            continue;
        end

        nremFrac = mean(nremLike(a:b));

        if nremFrac<CFG.min_nrem_like_fraction
            continue;
        end

        bpSeg = eegSp(a:b);

        minPeakDistance = max(1,round(fs/20));

        [~,locs] = findpeaks( ...
            bpSeg, ...
            'MinPeakDistance',minPeakDistance);

        nCycles = numel(locs);

        if nCycles<CFG.min_cycles || nCycles>CFG.max_cycles
            continue;
        end

        rawSeg = eeg(a:b);

        nfft = max(512,2^nextpow2(max(numel(rawSeg)*4,round(2*fs))));

        [Pxx,F] = periodogram( ...
            rawSeg, ...
            hann(numel(rawSeg)), ...
            nfft, ...
            fs);

        pSpin = peak_band(Pxx,F,CFG.spindle_band_hz);
        pLow = peak_band(Pxx,F,CFG.low_flank_hz);
        pHigh = peak_band(Pxx,F,CFG.high_flank_hz);

        if ~isfinite(pSpin) || ...
                ~isfinite(pLow) || ...
                ~isfinite(pHigh) || ...
                pSpin<=pLow || ...
                pSpin<=pHigh
            continue;
        end

        spinMask = F>=CFG.spindle_band_hz(1) & ...
                   F<=CFG.spindle_band_hz(2);

        Fs = F(spinMask);
        Ps = Pxx(spinMask);

        [~,q] = max(Ps);
        centralHz = Fs(q);

        [~,pkLocal] = max(abs(bpSeg));
        peakIndex = a+pkLocal-1;

        [~,energyPeakLocal] = max(eSeg);

        symmetry = NaN;
        if numel(eSeg)>1
            symmetry = (energyPeakLocal-1)/(numel(eSeg)-1);
        end

        tStart = (a-1)/fs;
        tEnd = b/fs;
        tPeak = (peakIndex-1)/fs;

        peakToPeak = max(bpSeg)-min(bpSeg);

        sessionPhase = condition_phase( ...
            condition,tPeak,mouseID,WARM,DRUG);

        warmingMatchedPhase = matched_phase( ...
            "warming",tPeak,mouseID,WARM,DRUG);

        drugMatchedPhase = matched_phase( ...
            "donepezil",tPeak,mouseID,WARM,DRUG);

        E = table( ...
            sid,mouseID,genotype,sex,cohort,cohortDisplay,condition, ...
            tStart,tEnd,tPeak,dur,centralHz,nCycles,symmetry,peakToPeak, ...
            nremFrac,max(eSeg),thrL,thrH,thrM, ...
            string(sessionPhase), ...
            string(warmingMatchedPhase), ...
            string(drugMatchedPhase), ...
            'VariableNames',{ ...
            'session_id','mouse_id','genotype','sex','cohort', ...
            'cohort_display','condition', ...
            'start_s','end_s','peak_s','duration_s', ...
            'central_frequency_hz','n_cycles','symmetry', ...
            'peak_to_peak','nrem_like_fraction','max_wavelet_energy', ...
            'threshold_low','threshold_high','threshold_max', ...
            'session_phase','warming_matched_phase','drug_matched_phase'});

        sessionEvents=[sessionEvents;E]; %#ok<AGROW>
    end

    EVENTS=[EVENTS;sessionEvents]; %#ok<AGROW>

    %% Window summaries

    AW = analysis_windows_for_session( ...
        condition,mouseID,durationS,WARM,DRUG);

    for w = 1:height(AW)

        t0 = max(0,double(AW.start_s(w)));
        t1 = min(durationS,double(AW.end_s(w)));

        if ~(isfinite(t0) && isfinite(t1) && t1>t0)
            continue;
        end

        i0 = max(1,floor(t0*fs)+1);
        i1 = min(N,ceil(t1*fs));

        nremS = sum(nremLike(i0:i1))/fs;

        ev = sessionEvents( ...
            sessionEvents.peak_s>=t0 & ...
            sessionEvents.peak_s<t1,:);

        nSp = height(ev);

        density = NaN;
        if nremS>0
            density=nSp/(nremS/60);
        end

        meanDur = NaN;
        medianDur = NaN;
        meanHz = NaN;
        meanCycles = NaN;
        meanAmp = NaN;

        if nSp>0
            meanDur=mean(ev.duration_s,'omitnan');
            medianDur=median(ev.duration_s,'omitnan');
            meanHz=mean(ev.central_frequency_hz,'omitnan');
            meanCycles=mean(ev.n_cycles,'omitnan');
            meanAmp=mean(ev.peak_to_peak,'omitnan');
        end

        W = table( ...
            sid,mouseID,genotype,sex,cohort,cohortDisplay,condition, ...
            string(AW.window(w)),t0,t1,nremS,nSp,density, ...
            meanDur,medianDur,meanHz,meanCycles,meanAmp, ...
            'VariableNames',{ ...
            'session_id','mouse_id','genotype','sex','cohort', ...
            'cohort_display','condition','window','start_s','end_s', ...
            'nrem_like_s','n_spindles','spindles_per_min_nrem', ...
            'mean_duration_s','median_duration_s', ...
            'mean_central_frequency_hz','mean_cycles', ...
            'mean_peak_to_peak'});

        WINDOWS=[WINDOWS;W]; %#ok<AGROW>
    end

    %% NREM transition analysis

    TR = build_transition_tables( ...
        scoreTime,scores,sessionEvents, ...
        sid,mouseID,genotype,sex,cohort,cohortDisplay,condition, ...
        WARM,DRUG,CFG);

    TRANSITIONS=[TRANSITIONS;TR.bouts]; %#ok<AGROW>
    TRANSITION_TC=[TRANSITION_TC;TR.timecourse]; %#ok<AGROW>

    %% QC row

    Q = table( ...
        sid,mouseID,genotype,sex,cohort,condition, ...
        string(rawRel),rawFs,fs,durationS,nremLikeS,height(sessionEvents), ...
        muE,sdE,thrL,thrH,thrM, ...
        'VariableNames',{ ...
        'session_id','mouse_id','genotype','sex','cohort','condition', ...
        'recording_file','raw_fs_hz','analysis_fs_hz','duration_s', ...
        'nrem_like_s','n_spindles','wavelet_energy_mean', ...
        'wavelet_energy_sd','threshold_low','threshold_high', ...
        'threshold_max'});

    QC=[QC;Q]; %#ok<AGROW>

    %% Limited deterministic QC plots

    if CFG.save_qc_plots && qcPlotCount<CFG.max_qc_plots

        step=max(1,floor(height(Sessions)/CFG.max_qc_plots));

        if mod(i-1,step)==0

            try
                make_qc_plot( ...
                    eeg,fs,scoreTime,scores,energy,sessionEvents, ...
                    thrL,thrH,sid,figDir);

                qcPlotCount=qcPlotCount+1;

            catch ME
                warning('QC figure failed for %s: %s',sid,ME.message);
            end
        end
    end

    fprintf('%4d spindles | %.2f/min NREM\n', ...
        height(sessionEvents), ...
        height(sessionEvents)/max(nremLikeS/60,eps));

end

%% ------------------------------------------------------------------------
% Mouse-level transition summaries
% -------------------------------------------------------------------------

TRANS_MOUSE = summarize_transition_bouts(TRANSITIONS);
TC_MOUSE = summarize_transition_timecourse(TRANSITION_TC);

%% ------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

writetable(EVENTS, ...
    fullfile(outDir,'spindle_events.csv'));

writetable(WINDOWS, ...
    fullfile(outDir,'spindle_window_summary.csv'));

writetable(TRANSITIONS, ...
    fullfile(outDir,'spindle_transition_bouts.csv'));

writetable(TRANS_MOUSE, ...
    fullfile(outDir,'spindle_transition_mouse_summary.csv'));

writetable(TRANSITION_TC, ...
    fullfile(outDir,'spindle_transition_timecourse_bouts.csv'));

writetable(TC_MOUSE, ...
    fullfile(outDir,'spindle_transition_timecourse_mouse.csv'));

writetable(QC, ...
    fullfile(outDir,'spindle_session_qc.csv'));

fprintf('\n============================================================\n');
fprintf('SPINDLE ANALYSIS COMPLETE\n');
fprintf('============================================================\n');
fprintf('Sessions analyzed: %d\n',height(QC));
fprintf('Detected spindles: %d\n',height(EVENTS));
fprintf('NREM transitions:  %d\n',height(TRANSITIONS));
fprintf('\nOutput:\n%s\n',outDir);
fprintf('============================================================\n\n');

%% ========================================================================
% Local functions
% =========================================================================

function tf = to_logical(x)

if islogical(x)
    tf=x;
elseif isnumeric(x)
    tf=x~=0;
else
    tf=ismember(lower(strtrim(string(x))), ...
        ["1","true","yes","y"]);
end
end


function score = reclassify_short_wake(score,maxWakeS)

score=double(score(:));

runs=state_runs(score);

for i=1:height(runs)

    if runs.state(i)==0

        n=runs.end_idx(i)-runs.start_idx(i)+1;

        if n<=maxWakeS
            score(runs.start_idx(i):runs.end_idx(i))=15;
        end
    end
end
end


function [stateSamples,validSamples] = scores_to_samples(timeS,scores,N,fs)

stateSamples=nan(N,1);
validSamples=false(N,1);

sampleTime=(0:N-1)'/fs;

% Each score is a 1-s epoch beginning at timeS(k).
idx=floor(sampleTime)+1;

valid = idx>=1 & ...
        idx<=numel(scores) & ...
        sampleTime < (timeS(end)+1);

stateSamples(valid)=scores(idx(valid));
validSamples(valid)=true;
end


function energy = fbsp_energy_chunked( ...
    x,fs,wname,freqs,weightByFreq,chunkS,marginS)

assert(exist('cwt','file')==2, ...
    'Wavelet Toolbox cwt() is required.');

assert(exist('centfrq','file')==2, ...
    'Wavelet Toolbox centfrq() is required.');

cf=centfrq(wname);
scales=cf*fs./freqs;

N=numel(x);
energy=nan(N,1);

coreN=max(1,round(chunkS*fs));
marginN=max(1,round(marginS*fs));

for core0=1:coreN:N

    core1=min(N,core0+coreN-1);

    ext0=max(1,core0-marginN);
    ext1=min(N,core1+marginN);

    seg=x(ext0:ext1).';

    C=cwt(seg,scales,wname);

    P=abs(C).^2;

    if weightByFreq
        P=P.*freqs(:);
    end

    e=mean(P,1);

    q0=core0-ext0+1;
    q1=q0+(core1-core0);

    energy(core0:core1)=e(q0:q1).';
end
end


function runs = logical_runs(mask)

mask=logical(mask(:));

d=diff([false;mask;false]);

s=find(d==1);
e=find(d==-1)-1;

runs=[s e];
end


function out = merge_intervals(in)

if isempty(in)
    out=zeros(0,2);
    return;
end

in=sortrows(in,1);
out=in(1,:);

for i=2:size(in,1)

    if in(i,1)<=out(end,2)+1
        out(end,2)=max(out(end,2),in(i,2));
    else
        out(end+1,:)=in(i,:); %#ok<AGROW>
    end
end
end


function p = peak_band(P,F,band)

idx=F>=band(1) & F<=band(2);

if any(idx)
    p=max(P(idx));
else
    p=NaN;
end
end


function rawRel = infer_raw_recording_path(sessionRow,dataRoot)

rawRel="";

vars=string(sessionRow.Properties.VariableNames);

if ~ismember("score_file",vars)
    return;
end

scoreRel=string(sessionRow.score_file(1));

if ismissing(scoreRel) || strlength(strtrim(scoreRel))==0
    return;
end

[folder,name,~]=fileparts(char(scoreRel));

base=regexprep(string(name),'_final_scoring_dash$','');

if endsWith(base,".mat")
    rawName=base;
else
    rawName=base+".mat";
end

candidate=string(fullfile(folder,char(rawName)));

if isfile(fullfile(dataRoot,char(candidate)))
    rawRel=candidate;
    return;
end

% Known June acquisition filename typo.
sid="";
if ismember("session_id",vars)
    sid=string(sessionRow.session_id(1));
end

if sid=="20260605_mouse63_warming_2"
    candidate="scored_files_june26batch/050626-m63-ammbt-2.mat";

    if isfile(fullfile(dataRoot,char(candidate)))
        rawRel=candidate;
        return;
    end
end
end


function W = analysis_windows_for_session( ...
    condition,mouseID,durationS,WARM,DRUG)

W=table( ...
    strings(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames',{'window','start_s','end_s'});

if condition=="baseline"
    W=[W; {"BASELINE_ZT3_9",0,min(durationS,6*3600)}];
end

[wStart,wEnd]=warming_bounds(mouseID,WARM);

if isfinite(wStart)

    if ~isfinite(wEnd) || wEnd<=wStart
        wEnd=wStart+3*3600;
    end

    pre0=max(0,wStart-3*3600);

    if condition=="warming"
        W=[W; {"WARM_PRE",pre0,wStart}];
        W=[W; {"WARM_ACTIVE",wStart,wEnd}];

    elseif condition=="baseline"
        W=[W; {"WARM_MATCH_PRE",pre0,wStart}];
        W=[W; {"WARM_MATCH_ACTIVE",wStart,wEnd}];
    end
end

d=drug_time(mouseID,DRUG);

if isfinite(d)

    pre0=max(0,d-3*3600);

    if condition=="donepezil"
        W=[W; {"DRUG_PRE",pre0,d}];
        W=[W; {"DRUG_0_3H",d,d+3*3600}];

    elseif condition=="baseline"
        W=[W; {"DRUG_MATCH_PRE",pre0,d}];
        W=[W; {"DRUG_MATCH_0_3H",d,d+3*3600}];
    end
end
end


function phase = condition_phase(condition,t,mouseID,WARM,DRUG)

phase="OTHER";

if condition=="baseline"

    if t>=0 && t<6*3600
        phase="BASELINE_ZT3_9";
    end

elseif condition=="warming"

    [s,e]=warming_bounds(mouseID,WARM);

    if isfinite(s)

        if ~isfinite(e) || e<=s
            e=s+3*3600;
        end

        if t>=max(0,s-3*3600) && t<s
            phase="WARM_PRE";
        elseif t>=s && t<e
            phase="WARM_ACTIVE";
        end
    end

elseif condition=="donepezil"

    d=drug_time(mouseID,DRUG);

    if isfinite(d)

        if t>=max(0,d-3*3600) && t<d
            phase="DRUG_PRE";
        elseif t>=d && t<d+3*3600
            phase="DRUG_0_3H";
        end
    end
end
end


function phase = matched_phase(kind,t,mouseID,WARM,DRUG)

phase="OTHER";

if kind=="warming"

    [s,e]=warming_bounds(mouseID,WARM);

    if ~isfinite(s)
        return;
    end

    if ~isfinite(e) || e<=s
        e=s+3*3600;
    end

    if t>=max(0,s-3*3600) && t<s
        phase="PRE";
    elseif t>=s && t<e
        phase="ACTIVE";
    end

elseif kind=="donepezil"

    d=drug_time(mouseID,DRUG);

    if ~isfinite(d)
        return;
    end

    if t>=max(0,d-3*3600) && t<d
        phase="PRE";
    elseif t>=d && t<d+3*3600
        phase="ACTIVE";
    end
end
end


function [s,e] = warming_bounds(mouseID,T)

s=NaN;
e=NaN;

if isempty(T)
    return;
end

if ~ismember("mouse_id",string(T.Properties.VariableNames))
    return;
end

Q=T(string(T.mouse_id)==mouseID,:);

if isempty(Q)
    return;
end

if ismember("include_analysis",string(Q.Properties.VariableNames))
    keep=to_logical(Q.include_analysis);
    Q=Q(keep,:);
end

if isempty(Q)
    return;
end

s=first_numeric_column(Q, ...
    ["warming_start_elapsed_s","warm_start_elapsed_s"], ...
    ["warming","start"], ...
    ["hms","source"]);

e=first_numeric_column(Q, ...
    ["warming_end_elapsed_s","warm_end_elapsed_s"], ...
    ["warming","end"], ...
    ["hms","source"]);
end


function d = drug_time(mouseID,T)

d=NaN;

if isempty(T)
    return;
end

if ~ismember("mouse_id",string(T.Properties.VariableNames))
    return;
end

Q=T(string(T.mouse_id)==mouseID,:);

if isempty(Q)
    return;
end

if ismember("include_analysis",string(Q.Properties.VariableNames))
    keep=to_logical(Q.include_analysis);
    Q=Q(keep,:);
end

if isempty(Q)
    return;
end

d=first_numeric_column(Q, ...
    ["drug_elapsed_s","drug_time_elapsed_s", ...
     "donepezil_elapsed_s","donepezil_time_elapsed_s"], ...
    ["drug"], ...
    ["saline","hms","source"]);

if ~isfinite(d)
    d=first_numeric_column(Q, ...
        strings(0,1), ...
        ["donepezil"], ...
        ["saline","hms","source"]);
end
end


function value = first_numeric_column(T,exactNames,containsTokens,excludeTokens)

value=NaN;

vars=string(T.Properties.VariableNames);
low=lower(vars);

for k=1:numel(exactNames)

    ix=find(low==lower(exactNames(k)),1);

    if ~isempty(ix)
        value=scalar_number(T.(vars(ix))(1));

        if isfinite(value)
            return;
        end
    end
end

candidate=true(size(vars));

for k=1:numel(containsTokens)
    candidate=candidate & contains(low,lower(containsTokens(k)));
end

for k=1:numel(excludeTokens)
    candidate=candidate & ~contains(low,lower(excludeTokens(k)));
end

ix=find(candidate);

for k=1:numel(ix)

    value=scalar_number(T.(vars(ix(k)))(1));

    if isfinite(value)
        return;
    end
end
end


function x = scalar_number(v)

x=NaN;

if isnumeric(v) || islogical(v)

    if ~isempty(v)
        x=double(v(1));
    end

else

    s=string(v);

    if ~ismissing(s) && strlength(strtrim(s))>0
        x=str2double(s);
    end
end
end


function OUT = build_transition_tables( ...
    scoreTime,scores,EV,sid,mouseID,genotype,sex,cohort,cohortDisplay, ...
    condition,WARM,DRUG,CFG)

OUT.bouts=table();
OUT.timecourse=table();

% Assimilate MA into NREM.
macro=scores(:);
macro(macro==15)=1;

runs=state_runs(macro);

boutCounter=0;

for i=1:height(runs)

    if runs.state(i)~=1
        continue;
    end

    e=runs.end_idx(i);

    if e>=numel(macro)
        continue;
    end

    nextState=macro(e+1);

    if nextState==2
        transitionType="NREM_to_REM";

    elseif nextState==0
        transitionType="NREM_to_WAKE";

    else
        continue;
    end

    s=runs.start_idx(i);

    tStart=scoreTime(s);
    tTransition=scoreTime(e+1);

    nremBoutDuration=tTransition-tStart;

    if ~(isfinite(nremBoutDuration) && nremBoutDuration>0)
        continue;
    end

    boutCounter=boutCounter+1;

    % Primary 25-s endpoint. If the NREM bout is shorter than 25 s, the
    % denominator is the available NREM duration and is exported explicitly.
    t25=max(tStart,tTransition-CFG.transition_primary_s);

    ev25=EV( ...
        EV.peak_s>=t25 & ...
        EV.peak_s<tTransition,:);

    effective25=tTransition-t25;

    rate25=height(ev25)/max(effective25/60,eps);

    % Last spindle within descriptive final 60 s.
    t60=max(tStart,tTransition-CFG.transition_timecourse_s);

    ev60=EV( ...
        EV.peak_s>=t60 & ...
        EV.peak_s<tTransition,:);

    lastPeakDistance=NaN;

    if ~isempty(ev60)
        lastPeakDistance=tTransition-max(ev60.peak_s);
    end

    meanDur=NaN;
    meanHz=NaN;
    meanCycles=NaN;

    if ~isempty(ev25)
        meanDur=mean(ev25.duration_s,'omitnan');
        meanHz=mean(ev25.central_frequency_hz,'omitnan');
        meanCycles=mean(ev25.n_cycles,'omitnan');
    end

    sessionPhase=condition_phase( ...
        condition,tTransition,mouseID,WARM,DRUG);

    warmingMatchedPhase=matched_phase( ...
        "warming",tTransition,mouseID,WARM,DRUG);

    drugMatchedPhase=matched_phase( ...
        "donepezil",tTransition,mouseID,WARM,DRUG);

    B=table( ...
        sid,mouseID,genotype,sex,cohort,cohortDisplay,condition, ...
        boutCounter,transitionType,tStart,tTransition,nremBoutDuration, ...
        effective25,height(ev25),rate25,lastPeakDistance, ...
        meanDur,meanHz,meanCycles, ...
        string(sessionPhase), ...
        string(warmingMatchedPhase), ...
        string(drugMatchedPhase), ...
        'VariableNames',{ ...
        'session_id','mouse_id','genotype','sex','cohort', ...
        'cohort_display','condition','bout_index','transition_type', ...
        'nrem_bout_start_s','transition_s','nrem_bout_duration_s', ...
        'primary_window_s','n_spindles_last25', ...
        'spindle_rate_last25_per_min', ...
        'last_spindle_peak_distance_s', ...
        'mean_spindle_duration_last25_s', ...
        'mean_spindle_frequency_last25_hz', ...
        'mean_spindle_cycles_last25', ...
        'session_phase','warming_matched_phase','drug_matched_phase'});

    OUT.bouts=[OUT.bouts;B]; %#ok<AGROW>

    % Descriptive 1-Hz occupancy in final 60 s.
    for rel=(-CFG.transition_timecourse_s):-1

        sec0=tTransition+rel;
        sec1=sec0+1;

        occupied=NaN;

        if sec0>=tStart

            q=EV.start_s<sec1 & ...
              EV.end_s>sec0;

            occupied=double(any(q));
        end

        T=table( ...
            sid,mouseID,genotype,sex,cohort,cohortDisplay,condition, ...
            boutCounter,transitionType,tTransition,rel,occupied, ...
            string(sessionPhase), ...
            string(warmingMatchedPhase), ...
            string(drugMatchedPhase), ...
            'VariableNames',{ ...
            'session_id','mouse_id','genotype','sex','cohort', ...
            'cohort_display','condition','bout_index','transition_type', ...
            'transition_s','relative_second','spindle_occupied', ...
            'session_phase','warming_matched_phase','drug_matched_phase'});

        OUT.timecourse=[OUT.timecourse;T]; %#ok<AGROW>
    end
end
end


function R = state_runs(x)

x=x(:);

change=[true;diff(x)~=0];

s=find(change);
e=[s(2:end)-1;numel(x)];

R=table( ...
    x(s),s,e, ...
    'VariableNames',{'state','start_idx','end_idx'});
end


function T = summarize_transition_bouts(B)

T=table();

if isempty(B)
    return;
end

keys=[ ...
    "mouse_id","genotype","sex","cohort","cohort_display","condition", ...
    "transition_type","session_phase", ...
    "warming_matched_phase","drug_matched_phase"];

[G,K]=findgroups(B(:,keys));

nB=splitapply(@numel,B.bout_index,G);

meanRate=splitapply( ...
    @(x) mean(x,'omitnan'), ...
    B.spindle_rate_last25_per_min,G);

medianRate=splitapply( ...
    @(x) median(x,'omitnan'), ...
    B.spindle_rate_last25_per_min,G);

meanLast=splitapply( ...
    @(x) mean(x,'omitnan'), ...
    B.last_spindle_peak_distance_s,G);

meanDur=splitapply( ...
    @(x) mean(x,'omitnan'), ...
    B.mean_spindle_duration_last25_s,G);

meanHz=splitapply( ...
    @(x) mean(x,'omitnan'), ...
    B.mean_spindle_frequency_last25_hz,G);

meanCycles=splitapply( ...
    @(x) mean(x,'omitnan'), ...
    B.mean_spindle_cycles_last25,G);

T=K;

T.n_transition_bouts=nB;
T.mean_spindle_rate_last25_per_min=meanRate;
T.median_spindle_rate_last25_per_min=medianRate;
T.mean_last_spindle_peak_distance_s=meanLast;
T.mean_spindle_duration_last25_s=meanDur;
T.mean_spindle_frequency_last25_hz=meanHz;
T.mean_spindle_cycles_last25=meanCycles;
end


function T = summarize_transition_timecourse(B)

T=table();

if isempty(B)
    return;
end

keys=[ ...
    "mouse_id","genotype","sex","cohort","cohort_display","condition", ...
    "transition_type","session_phase", ...
    "warming_matched_phase","drug_matched_phase","relative_second"];

[G,K]=findgroups(B(:,keys));

meanOcc=splitapply( ...
    @(x) mean(x,'omitnan'), ...
    B.spindle_occupied,G);

nB=splitapply( ...
    @(x) sum(isfinite(x)), ...
    B.spindle_occupied,G);

T=K;
T.mean_spindle_occupancy=meanOcc;
T.n_bouts_contributing=nB;
end


function make_qc_plot( ...
    eeg,fs,scoreTime,scores,energy,EV,thrL,thrH,sid,figDir)

N=numel(eeg);
winN=min(N,round(120*fs));

[state,valid]=scores_to_samples(scoreTime,scores,N,fs);

nrem=valid & (state==1 | state==15);

if N>winN
    richness=conv(double(nrem),ones(winN,1),'valid');
    [~,i0]=max(richness);
else
    i0=1;
end

i1=min(N,i0+winN-1);

t=((i0:i1)-1)'/fs;

fig=figure( ...
    'Color','w', ...
    'Visible','off', ...
    'Position',[80 80 1200 650]);

tl=tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax=nexttile(tl);
plot(ax,t,eeg(i0:i1),'k');
ylabel(ax,'EEG');
title(ax,strrep(char(sid),'_','\_'));
hold(ax,'on');

q=EV.end_s>=t(1) & EV.start_s<=t(end);

for j=find(q)'
    xline(ax,EV.start_s(j),'-');
    xline(ax,EV.end_s(j),':');
end

ax=nexttile(tl);
plot(ax,t,energy(i0:i1),'k');
hold(ax,'on');

yline(ax,thrL,'--');
yline(ax,thrH,'-');

ylabel(ax,'Fbsp energy');
xlabel(ax,'Time (s)');

exportgraphics( ...
    fig, ...
    fullfile(figDir,char(sid+"_spindle_qc.png")), ...
    'Resolution',180);

close(fig);
end

function y = apply_notch_compat(x,fs,f0,Q)
%APPLY_NOTCH_COMPAT Notch filter compatible with older MATLAB versions.
%
% Uses iirnotch when available. Otherwise uses a narrow Butterworth
% band-stop filter with approximately the same bandwidth (f0/Q).

x = double(x(:));

if ~(isfinite(f0) && isfinite(Q) && f0>0 && Q>0 && f0<fs/2)
    y = x;
    return;
end

if exist('iirnotch','file')==2

    wo = f0/(fs/2);
    bw = wo/Q;

    [b,a] = iirnotch(wo,bw);

else

    bwHz = f0/Q;

    f1 = max(0.01, f0-bwHz/2);
    f2 = min(fs/2-0.01, f0+bwHz/2);

    if f1>=f2
        y = x;
        return;
    end

    [b,a] = butter(2,[f1 f2]/(fs/2),'stop');

end

y = filtfilt(b,a,x);

end

