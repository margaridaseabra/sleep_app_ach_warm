%% RUN_14_SPINDLE_PROBABILITY
% Reproduce/extend thesis-style pre-transition spindle occupancy analysis.
% Requires existing outputs of run_10_spindle_analysis.m; NO EEG redetection.
%
% Expected CSVs under P.eeg_output/spindles:
%   spindle_transition_bouts.csv
%   spindle_transition_timecourse_bouts.csv
%
% Produces under outputs/eeg/spindles/probability:
%   qc_summary.txt
%   eligible_bouts_25s.csv
%   mouse_probability_summary_25s.csv
%   mouse_probability_timecourse_60s.csv
%   figures/spindle_probability_{all,thesis,new}.{png,pdf,fig}
%
% Scientific conventions:
% - Each occupancy value is 1 if any part of a spindle overlaps that SECOND.
% - -25:-16 is the 10 one-second bins spanning [-25,-15) seconds.
% - -25:-1 is the 25 one-second bins spanning [-25,0) seconds.
% - A bout enters a scalar endpoint only when all 25 bins are present.
% - A bout enters 60-s curves only when all 60 bins are present.
% - The complete history window must fall inside the specified phase.
% - Mouse (not bout) is the group-level replication unit.
% - Session clock is assumed to start at ZT3: t=0 at ZT3.
%
% Compatible with MATLAB table/string operations (R2021a+ recommended).

clear; clc;
scriptFile = mfilename('fullpath');
if isempty(scriptFile), error('Save this as run_14_spindle_probability.m and run the saved script.'); end
repoDir = fileparts(fileparts(scriptFile));
addpath(fullfile(repoDir,'config'));
if exist('paths_config','file') ~= 2
    error('Cannot find config/paths_config.m. Put this script in your repository pipelines/ folder.');
end
P = paths_config();
inDir = fullfile(P.eeg_output,'spindles');
outDir = fullfile(inDir,'probability');
figDir = fullfile(outDir,'figures');
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(figDir), mkdir(figDir); end

bFile = fullfile(inDir,'spindle_transition_bouts.csv');
tcFile = fullfile(inDir,'spindle_transition_timecourse_bouts.csv');
assert(isfile(bFile), 'Missing %s. Run the existing spindle detector first.', bFile);
assert(isfile(tcFile), 'Missing %s. Run the existing spindle detector first.', tcFile);

MIN_BOUTS_FOR_PLOTS = 3;
H = 3600; % Seconds after ZT3 at recording start.

fprintf('\n=== SPINDLE PROBABILITY / EXISTING EVENT CSVs ===\n');
B = readtable(bFile, 'TextType','string');
TC = readtable(tcFile, 'TextType','string');
require_vars(B, {'session_id','bout_index','mouse_id','genotype', ...
    'cohort','condition','transition_type','transition_s','nrem_bout_duration_s'}, bFile);
require_vars(TC, {'session_id','bout_index','relative_second', ...
    'spindle_occupied'}, tcFile);

% Join by SESSION + BOUT rather than bout_index alone (which repeats).
B.session_id = string(B.session_id);
B.bout_index = string(B.bout_index);
B.mouse_id = string(B.mouse_id);
B.genotype = upper(strtrim(string(B.genotype)));
B.genotype(B.genotype == "APP/PS1") = "APP";
B.cohort = lower(strtrim(string(B.cohort)));
B.condition = lower(strtrim(string(B.condition)));
B.transition_type = string(B.transition_type);
B.transition_s = double(B.transition_s);
B.nrem_bout_duration_s = double(B.nrem_bout_duration_s);

TC.session_id = string(TC.session_id);
TC.bout_index = string(TC.bout_index);
TC.relative_second = double(TC.relative_second);
TC.spindle_occupied = double(TC.spindle_occupied);

B.key = B.session_id + "||" + B.bout_index;
TC.key = TC.session_id + "||" + TC.bout_index;
if numel(unique(B.key)) ~= height(B)
    error('Duplicate (session_id, bout_index) keys in spindle_transition_bouts.csv.');
end

% Do not silently accept unknown state/condition/genotype names.
goodLabels = ismember(B.genotype,["WT","APP"]) & ...
    ismember(B.condition,["baseline","warming","donepezil"]) & ...
    ismember(B.transition_type,["NREM_to_REM","NREM_to_WAKE"]);
if nnz(~goodLabels)>0
    warning('Ignoring %d bouts with unexpected genotype, condition or transition labels.',nnz(~goodLabels));
    B = B(goodLabels,:);
end
B.source_idx = (1:height(B))';

% Build a fixed [bout x 60 second] matrix. NaN means unavailable, not no spindle.
M = nan(height(B),60);
TC = TC(isfinite(TC.relative_second) & ...
        abs(TC.relative_second-round(TC.relative_second))<1e-8 & ...
        TC.relative_second>=-60 & TC.relative_second<=-1,:);
[isKnown, boutRow] = ismember(TC.key,B.key);
TC = TC(isKnown,:);
boutRow = boutRow(isKnown);
secondColumn = TC.relative_second + 61; % -60 -> 1; -1 -> 60
linearIdx = sub2ind(size(M),boutRow,secondColumn);
if numel(unique(linearIdx)) ~= numel(linearIdx)
    error('Duplicate bout/relative_second pairs in timecourse CSV. Verify run_10 output.');
end
if any(~ismember(TC.spindle_occupied(isfinite(TC.spindle_occupied)),[0 1]))
    error('spindle_occupied must be binary 0/1 (or NaN for missing).');
end
M(linearIdx) = TC.spindle_occupied;

% Keep only complete NREM windows. The detector's primary_window_s is an
% additional QC condition when exported, but complete occupancy is essential.
validDuration25 = B.nrem_bout_duration_s >= 25-1e-6;
validDuration60 = B.nrem_bout_duration_s >= 60-1e-6;
if ismember('primary_window_s',B.Properties.VariableNames)
    validDuration25 = validDuration25 & double(B.primary_window_s)>=25-1e-6;
end
full25 = validDuration25 & all(isfinite(M(:,36:60)),2);
full60 = validDuration60 & all(isfinite(M),2);
fprintf('Source bouts: %d | full 25 s: %d | full 60 s: %d\n', ...
    height(B), nnz(full25), nnz(full60));

% Define exactly which intervention phase contains the HISTORY, not merely
% the scored transition boundary. Includes supplementary 3rd-hour blocks.
B.phase = repmat("",height(B),1);
B.phase_start_s = nan(height(B),1);
B.phase_end_s = nan(height(B),1);
for i=1:height(B)
    t = B.transition_s(i);
    if ~(isfinite(t) && t>=0 && t<9*H), continue; end
    if t < 3*H
        B.phase_start_s(i)=0; B.phase_end_s(i)=3*H;
        switch B.condition(i)
            case "baseline",  B.phase(i)="BASELINE_ZT3_6";
            case "warming",   B.phase(i)="PREWARM_ZT3_6";
            case "donepezil", B.phase(i)="PREDRUG_ZT3_6";
        end
    elseif t < 6*H
        B.phase_start_s(i)=3*H; B.phase_end_s(i)=6*H;
        switch B.condition(i)
            case "baseline",  B.phase(i)="BASELINE_ZT6_9";
            case "warming",   B.phase(i)="WARM_ACTIVE";
            case "donepezil", B.phase(i)="DRUG_0_3H";
        end
    else
        B.phase_start_s(i)=6*H; B.phase_end_s(i)=9*H;
        switch B.condition(i)
            case "baseline",  B.phase(i)="BASELINE_ZT9_12";
            case "warming",   B.phase(i)="WARM_COOLDOWN";
            case "donepezil", B.phase(i)="DRUG_3_6H";
        end
    end
end

% Add a pooled baseline ZT3-9 reference without treating this as an
% intervention phase boundary. This duplicates baseline observations into
% an ALTERNATIVE group (not the same statistical sample counted twice).
pooled = B(B.condition=="baseline" & B.transition_s>=0 & ...
           B.transition_s<6*H,:);
pooled.phase(:)="BASELINE_ZT3_9";
pooled.phase_start_s(:)=0;
pooled.phase_end_s(:)=6*H;
S = [B; pooled];

% Primary per-bout endpoints; histories cannot cross condition boundaries.
idx25 = S.source_idx;
eligible25 = full25(idx25) & strlength(S.phase)>0 & ...
    S.transition_s >= S.phase_start_s+25 & ...
    S.transition_s < S.phase_end_s;
E = S(eligible25,:);
X25 = M(E.source_idx,36:60);
E.occupancy_probability_25_15 = mean(X25(:,1:10),2); % [-25,-15)
E.any_spindle_25 = double(any(X25>0,2));
E.occupancy_fraction_25 = mean(X25,2);
E = removevars(E,intersect({'key','source_idx','phase_start_s', ...
    'phase_end_s'},E.Properties.VariableNames));
writetable(E,fullfile(outDir,'eligible_bouts_25s.csv'));

% Group by mouse, NOT by bout, to avoid pseudoreplication.
keys = {'mouse_id','genotype','cohort','condition','phase','transition_type'};
[G,K] = findgroups(E(:,keys));
Mouse = K;
Mouse.n_bouts = splitapply(@numel,E.any_spindle_25,G);
Mouse.occupancy_probability_25_15 = splitapply(@mean,E.occupancy_probability_25_15,G);
Mouse.any_spindle_25 = splitapply(@mean,E.any_spindle_25,G);
Mouse.occupancy_fraction_25 = splitapply(@mean,E.occupancy_fraction_25,G);
writetable(Mouse,fullfile(outDir,'mouse_probability_summary_25s.csv'));

% QC: How many eligible bouts and animals contribute to each group?
qcKeys={'cohort','genotype','condition','phase','transition_type'};
[Gqc,Kqc]=findgroups(E(:,qcKeys));
QcGroups=Kqc;
QcGroups.n_bouts=splitapply(@numel,E.any_spindle_25,Gqc);
[Gqm,Kqm]=findgroups(Mouse(:,qcKeys));
QcMice=Kqm;
QcMice.n_mice_any=splitapply(@numel,Mouse.n_bouts,Gqm);
QcMice.n_mice_3plus=splitapply(@(x) nnz(x>=3),Mouse.n_bouts,Gqm);
QcGroups=outerjoin(QcGroups,QcMice,'Keys',qcKeys,'MergeKeys',true);
writetable(QcGroups,fullfile(outDir,'qc_counts_by_group.csv'));

% Only full, non-straddling 60-s histories enter the timecourse figure.
idx60 = S.source_idx;
eligible60 = full60(idx60) & strlength(S.phase)>0 & ...
    S.transition_s >= S.phase_start_s+60 & ...
    S.transition_s < S.phase_end_s;
R = S(eligible60,:);
[G60,K60] = findgroups(R(:,keys));
nGroups=height(K60);
Y=nan(nGroups,60); nB=zeros(nGroups,1);
for k=1:nGroups
    q = G60==k;
    nB(k)=nnz(q);
    Y(k,:) = mean(M(R.source_idx(q),:),1);
end
Timecourse = K60(repelem((1:nGroups)',60),:);
Timecourse.relative_second = repmat((-60:-1)',nGroups,1);
Timecourse.spindle_probability = reshape(Y.',[],1);
Timecourse.n_bouts = repelem(nB,60);
writetable(Timecourse,fullfile(outDir,'mouse_probability_timecourse_60s.csv'));

% QC summary (retained counts can differ between endpoint types).
qcFile=fullfile(outDir,'qc_summary.txt');
fid=fopen(qcFile,'w');
fprintf(fid,'Spindle probability pipeline QC\n');
fprintf(fid,'Source unique bouts: %d\n',height(B));
fprintf(fid,'Full 25-s histories before phase exclusion: %d\n',nnz(full25));
fprintf(fid,'Full 60-s histories before phase exclusion: %d\n',nnz(full60));
fprintf(fid,'Eligible bout/phase records (includes alternative pooled BL): %d\n',height(E));
fprintf(fid,'Eligible 60-s bout/phase records: %d\n',height(R));
fprintf(fid,'Mouse-by-phase-by-transition scalar rows: %d\n',height(Mouse));
fprintf(fid,'Mouse-by-phase-by-transition full 60-s curves: %d\n',nGroups);
fprintf(fid,'Mouse-level summary rows with >=%d bouts: %d\n', ...
    MIN_BOUTS_FOR_PLOTS,nnz(Mouse.n_bouts>=MIN_BOUTS_FOR_PLOTS));
fprintf(fid,'Assumption: recording t=0 is ZT3 and phase boundary is at t=10800 s.\n');
fprintf(fid,'No results are inferential until run_15 is completed.\n');
fclose(fid);

% Six thesis-style panels by transition outcome x recording condition.
conditions=["baseline","warming","donepezil"];
phases=["BASELINE_ZT3_9","WARM_ACTIVE","DRUG_0_3H"];
transitions=["NREM_to_REM","NREM_to_WAKE"];
cohorts=["all","thesis","new"];
colors=[0.40 0.40 0.40; 0.17 0.45 0.75];
for c=1:numel(cohorts)
    f=figure('Color','w','Position',[80 80 1350 680]);
    axs=gobjects(6,1); largest=0;
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    for tr=1:2
        for cc=1:3
            pos=(tr-1)*3+cc;
            axs(pos)=nexttile; hold on;
            hasCurve=false;
            D=Timecourse(Timecourse.condition==conditions(cc) & ...
                Timecourse.phase==phases(cc) & ...
                Timecourse.transition_type==transitions(tr) & ...
                Timecourse.n_bouts>=MIN_BOUTS_FOR_PLOTS,:);
            if cohorts(c)~="all", D=D(D.cohort==cohorts(c),:); end
            for gi=1:2
                genotype=["WT","APP"];
                A=D(D.genotype==genotype(gi),:);
                mice=unique(A.mouse_id);
                W=nan(numel(mice),60);
                for mi=1:numel(mice)
                    d=A(A.mouse_id==mice(mi),:);
                    [tSort,ord]=sort(d.relative_second);
                    if numel(tSort)==60 && isequal(tSort,(-60:-1)')
                        W(mi,:)=d.spindle_probability(ord)';
                    end
                end
                W=W(all(isfinite(W),2),:);
                if isempty(W), continue; end
                % 5-s smoothing for visualization only, across mouse curves.
                W=movmean(W,5,2);
                mu=mean(W,1)*100;
                if size(W,1)>1
                    se=std(W,0,1)/sqrt(size(W,1))*100;
                else
                    se=zeros(1,60);
                end
                xx=-60:-1;
                lower=max(0,mu-se); upper=min(100,mu+se);
                fill([xx fliplr(xx)],[lower fliplr(upper)],colors(gi,:), ...
                    'FaceAlpha',0.16,'EdgeColor','none','HandleVisibility','off');
                hasCurve=true;
                plot(xx,mu,'Color',colors(gi,:),'LineWidth',1.9, ...
                    'DisplayName',sprintf('%s (n=%d mice)',char(genotype(gi)),size(W,1)));
                largest=max(largest,max(upper));
            end
            xline(-25,':','HandleVisibility','off');
            xline(-15,':','HandleVisibility','off');
            xlim([-60 0]); xlabel('Seconds before transition');
            ylabel('1-s spindle occupancy (%)');
            if tr==1, outcome='NREM to REM'; else, outcome='NREM to Wake'; end
            title(sprintf('%s | %s',char(conditions(cc)),outcome),'Interpreter','none');
            if hasCurve, legend('Location','best'); end
            box off;
        end
    end
    % A shared readable scale: occupancy is typically far below 100%.
    top=min(100,max(8,ceil(largest*1.15/2)*2));
    for k=1:numel(axs), ylim(axs(k),[0 top]); end
    name="spindle_probability_"+cohorts(c);
    exportgraphics(f,fullfile(figDir,char(name+'.png')),'Resolution',300);
    exportgraphics(f,fullfile(figDir,char(name+'.pdf')),'ContentType','vector');
    savefig(f,fullfile(figDir,char(name+'.fig')));
    close(f);
end
fprintf('Scalar mouse summaries: %d rows; timecourse mouse groups: %d\n',height(Mouse),nGroups);
fprintf('Outputs: %s\n',outDir);
fprintf('Next: run(''pipelines/run_15_spindle_probability_stats.m'')\n');

function require_vars(T,names,file)
missing=setdiff(names,T.Properties.VariableNames);
if ~isempty(missing)
    error('CSV %s is missing columns: %s\nFound: %s',file, ...
        strjoin(missing,', '),strjoin(T.Properties.VariableNames,', '));
end
end
