%% RUN_22_PAC_REAL_COMODULOGRAMS
% Recalculate REAL frequency-resolved baseline Morlet PAC maps from raw EEG.
% Uses the exact accepted baseline windows / frozen NREM transitions from
% run_20_morlet_pac_baseline_full.m, with its locked bandwidth calibration,
% preprocessing, sampling support and deterministic surrogate shifts.
%
% Per-window surrogate-normalized PAC Z map --> mean over windows within a
% bout --> mean over bouts within a mouse --> mean over mice within genotype.
% Never pool all 10-s windows as independent mice.
%
% Maps are descriptive. DO NOT treat colored pixels as significant clusters.
% No changes to run_10, run_11*, run_20 or their outputs.
%
% Optional BEFORE running:
%   setenv('PAC_MAP_GRID','locked')  % same frequencies as run_20 (DEFAULT)
%   setenv('PAC_MAP_GRID','fine')    % more frequency centers, recomputed EEG
% Optional: setenv('PAC_MAP_FORCE','1') to regenerate per-session caches.
clear; clc;
repoDir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoDir,'config'));
P=paths_config();
C.fs=200; C.win_s=10; C.nsurr=200; C.min_shift_s=1;
C.seed=20261023; C.guard_s=3.2; C.padding_s=8;
C.edge_trim_s=1; C.wavelet_support_sigma=4; C.artifact_threshold=15;
mode=lower(string(getenv('PAC_MAP_GRID')));
if strlength(mode)==0,mode="locked";end
assert(any(mode==["locked","fine"]),'PAC_MAP_GRID must be locked or fine');
if mode=="locked"
    C.pf_nrem=1:4; C.af_nrem=9:16;
    C.pf_rem=4:2:10; C.af_rem=30:10:80;
else
    % Genuine newly evaluated frequencies -- NOT image interpolation.
    % Same Hilbert-filter FWHM calibration rule, continuous across centers.
    C.pf_nrem=1:0.25:4; C.af_nrem=9:0.5:16;
    C.pf_rem=4:0.5:10; C.af_rem=30:2:80;
end
C.cycles_pn=cycles_from_filter(C.pf_nrem,C.fs,0.5);
C.cycles_an=cycles_from_filter(C.af_nrem,C.fs,4);
C.cycles_pr=cycles_from_filter(C.pf_rem,C.fs,0.5);
C.cycles_ar=cycles_from_filter(C.af_rem,C.fs,10);
root=fullfile(P.eeg_output,'pac','baseline_morlet_bandwidth_v20');
refFile=fullfile(root,'paired_pac_windows.csv');
refMouseFile=fullfile(root,'paired_pac_mouse.csv');
transitionFile=fullfile(P.eeg_output,'spindles','spindle_transition_bouts.csv');
scoreDir=fullfile(P.sleep_output,'all','standardized_scores');
assert(isfile(refFile)&&isfile(refMouseFile), ...
    'Run full Morlet baseline v2.0 first; its window and mouse CSVs are required.');
assert(isfile(transitionFile)&&isfolder(scoreDir), ...
    'Missing frozen transitions or standardized score directory.');
output=fullfile(P.eeg_output,'pac',['real_comodulograms_v22_' char(mode)]);
cache=fullfile(output,'per_session');
figDir=fullfile(output,'figures');
if ~isfolder(cache),mkdir(cache);end
if ~isfolder(figDir),mkdir(figDir);end
force=strcmp(getenv('PAC_MAP_FORCE'),'1');
Ref=readtable(refFile,'TextType','string','VariableNamingRule','preserve');
B=readtable(transitionFile,'TextType','string','VariableNamingRule','preserve');
S=readtable(P.sessions,'TextType','string','VariableNamingRule','preserve');
M=readtable(P.mice,'TextType','string','VariableNamingRule','preserve');
S=S(string(S.condition)=="baseline" & double(S.include_sleep)==1,:);
assert(height(S)==16,'Expected 16 baseline sessions, found %d',height(S));
keyCols={'session_id','context','bout_index','window_start_s'};
assert(height(unique(Ref(:,keyCols)))==height(Ref),'Duplicated reference windows');
assert(all(ismember(["NREM_to_REM","NREM_to_WAKE","REM"], ...
    unique(string(Ref.context)))),'Missing a PAC context in reference windows');

% Record exact analysis settings, including frequency grids and cycles.
Settings=struct('mode',char(mode),'fs_hz',C.fs,'n_surrogates',C.nsurr,...
    'seed',C.seed,'guard_s',C.guard_s,'padding_s',C.padding_s, ...
    'edge_trim_s',C.edge_trim_s,'hierarchy', ...
    'window Z -> bout mean -> mouse mean -> cohort/genotype mean -> equal-cohort group mean', ...
    'note','Maps are descriptive, NOT pixelwise statistical inference');
save(fullfile(output,'analysis_settings.mat'),'Settings','C');

fprintf('\n=== REAL MORLET COMODULOGRAMS | %s | 16 baseline mice ===\n',char(mode));
fprintf('Locked to existing accepted windows and spindle transitions; %d surrogates.\n',C.nsurr);
fprintf('All frequency cells are calculated from raw EEG, not interpolated images.\n');
if mode=="fine",fprintf('FINE grid: this recalculates more wavelet frequencies per EEG window.\n');end
BoutLong=table();QC=table();
for i=1:height(S)
    sid=string(S.session_id(i));mid=string(S.mouse_id(i));
    cohort=string(S.cohort(i));
    mr=M(string(M.mouse_id)==mid,:);
    assert(height(mr)==1,'No unique metadata for %s',char(mid));
    genotype=string(mr.genotype(1));
    h=Ref(string(Ref.session_id)==sid,:);
    assert(~isempty(h),'Missing reference windows for %s',char(sid));
    cacheB=fullfile(cache,char(sid+'_bout_maps.csv'));
    cacheQ=fullfile(cache,char(sid+'_qc.csv'));
    if ~force && isfile(cacheB) && isfile(cacheQ)
        t=readtable(cacheB,'TextType','string','VariableNamingRule','preserve');
        q=readtable(cacheQ,'TextType','string','VariableNamingRule','preserve');
        BoutLong=[BoutLong;t];QC=[QC;q]; %#ok<AGROW>
        fprintf('[%02d/16] %s: loaded cache\n',i,char(sid));
        continue;
    end
    raw=fullfile(P.data_root,char(S.recording_file(i)));
    sf=fullfile(scoreDir,char(sid+'.csv'));
    assert(isfile(raw)&&isfile(sf),'Raw EEG / scores missing for %s',char(sid));
    R=load(raw,'eeg','eeg_frequency');
    eeg=double(R.eeg(:));rawFs=double(R.eeg_frequency);
    assert(isscalar(rawFs)&&isfinite(rawFs)&&rawFs>0,'Invalid sample rate');
    if abs(rawFs-C.fs)>1e-9
        [nn,dd]=rat(C.fs/rawFs,1e-10);
        eeg=resample(eeg,nn,dd);
    end
    eeg=eeg-mean(eeg,'omitnan');
    assert(all(isfinite(eeg)),'Non-finite EEG in %s',char(sid));
    if exist('iirnotch','file')==2
        [bb,aa]=iirnotch(50/(C.fs/2),(50/35)/(C.fs/2));
    else
        [bb,aa]=butter(2,[(50-50/70),(50+50/70)]/(C.fs/2),'stop');
    end
    eeg=filtfilt(bb,aa,eeg);
    Score=readtable(sf,'TextType','string','VariableNamingRule','preserve');
    ts=double(Score.time_s(:));state=double(Score.score(:));
    v=isfinite(ts)&isfinite(state);ts=ts(v);state=state(v);ts=ts-ts(1);
    state=reclassify_short_wake(state,15);
    state(state==15)=1;
    runs=state_runs(state);
    remRuns=runs(runs.state==2,:);
    br=B(string(B.session_id)==sid,:);
    perBout=table();nChecked=0;nSkipped=0;maxReferenceError=0;
    for ctx=["NREM_to_REM","NREM_to_WAKE","REM"]
        hh=h(string(h.context)==ctx,:);
        if isempty(hh),continue;end
        ids=sort(unique(double(hh.bout_index)));
        if ctx=="REM"
            pf=C.pf_rem;af=C.af_rem;
        else
            pf=C.pf_nrem;af=C.af_nrem;
        end
        [AMP,PHASE]=ndgrid(af,pf);nCells=numel(AMP);
        for jj=1:numel(ids)
            bid=ids(jj);
            w=hh(double(hh.bout_index)==bid,:);
            [~,ord]=sort(double(w.window_start_s));w=w(ord,:);
            if ctx~="REM"
                q=br(double(br.bout_index)==bid & string(br.transition_type)==ctx,:);
                assert(height(q)==1,'NREM transition ID mismatch for %s',char(sid));
                t1=double(q.transition_s(1));
                t0=max(double(q.nrem_bout_start_s(1)),t1-70);
                assert(all(double(w.window_start_s)>=t1-60-1e-3 & ...
                    double(w.window_start_s)+10<=t1+1e-3), ...
                    'Reference window outside pre-transition interval');
            else
                assert(bid>=1&&bid<=height(remRuns),'REM bout ID mismatch');
                t0=ts(remRuns.start_idx(bid));
                t1=ts(remRuns.end_idx(bid))+1;
                t0=max(t0,min(double(w.window_start_s))-10);
                t1=min(t1,max(double(w.window_start_s))+20);
            end
            a=max(1,round(t0*C.fs)+1);
            b=min(numel(eeg),round(t1*C.fs));
            if b<=a+C.fs*10
                nSkipped=nSkipped+height(w);continue;
            end
            xx=eeg(a:b).';
            tStart=(a-1)/C.fs;
            [phase,amp]=morlet_components(xx,C.fs,pf,af,ctx,C);
            time=tStart+(0:numel(xx)-1)/C.fs;
            sumMap=zeros(size(AMP));count=0;
            for j=1:height(w)
                ws=double(w.window_start_s(j));
                use=time>=ws & time<ws+C.win_s & ...
                    time>=tStart+C.guard_s & ...
                    time<time(end)-C.guard_s & ...
                    time>=ws+C.edge_trim_s & ...
                    time<ws+C.win_s-C.edge_trim_s;
                nSamples=sum(use);
                if nSamples<5*C.fs
                    nSkipped=nSkipped+1;continue;
                end
                % Identical sample support to v2.0, including truncated windows.
                expectedN=double(w.n_common_samples(j));
                assert(nSamples==expectedN, ...
                    'Sample support changed from v2.0: %s, %.3f s, %d vs %d', ...
                    char(sid),ws,nSamples,expectedN);
                rng(window_seed(sid,ctx,bid,ws,C.seed),'twister');
                shifts=randi([round(C.fs*C.min_shift_s), ...
                    nSamples-round(C.fs*C.min_shift_s)],C.nsurr,1);
                Z=surrogate_zmap(phase(:,use),amp(:,use),shifts,C);
                if mode=="locked"
                    errorZ=abs(mean(Z(:))-double(w.morlet_roi_mean_z(j)));
                    maxReferenceError=max(maxReferenceError,errorZ);
                end
                sumMap=sumMap+Z;
                count=count+1;nChecked=nChecked+1;
            end
            if count==0,continue;end
            zBout=sumMap/count;
            t=table(repmat(sid,nCells,1),repmat(mid,nCells,1), ...
                repmat(genotype,nCells,1),repmat(cohort,nCells,1), ...
                repmat(ctx,nCells,1),repmat(bid,nCells,1), ...
                repmat(count,nCells,1),PHASE(:),AMP(:),zBout(:), ...
                'VariableNames',{'session_id','mouse_id','genotype','cohort', ...
                'context','bout_index','n_windows','phase_hz', ...
                'amp_hz','pac_z'});
            perBout=[perBout;t]; %#ok<AGROW>
        end
    end
    if mode=="locked"
        assert(maxReferenceError<1e-3, ...
            'Map calculation differs from v2.0 scalar Z; max abs error %g', ...
            maxReferenceError);
    end
    q=table(sid,mid,genotype,cohort,height(h),nChecked,nSkipped, ...
        maxReferenceError,rawFs, ...
        'VariableNames',{'session_id','mouse_id','genotype','cohort', ...
        'reference_windows','evaluated_windows','skipped_windows', ...
        'max_absolute_reference_z_error','raw_fs_hz'});
    writetable(perBout,cacheB);writetable(q,cacheQ);
    BoutLong=[BoutLong;perBout];QC=[QC;q]; %#ok<AGROW>
    fprintf('[%02d/16] %-34s %d/%d matched maps, skipped %d\n', ...
        i,char(sid),nChecked,height(h),nSkipped);
    clear eeg R;
end
assert(~isempty(BoutLong),'No maps generated');
writetable(BoutLong,fullfile(output,'bout_comodulograms_long.csv'));
writetable(QC,fullfile(output,'pac_map_qc.csv'));
mouseKeys={'session_id','mouse_id','genotype','cohort','context','phase_hz','amp_hz'};
[gm,km]=findgroups(BoutLong(:,mouseKeys));
Mouse=km;
Mouse.pac_z=splitapply(@mean,BoutLong.pac_z,gm);
Mouse.n_bouts=splitapply(@numel,BoutLong.pac_z,gm);
writetable(Mouse,fullfile(output,'mouse_comodulograms_long.csv'));
% Cohort-balanced figures: average mice within genotype & cohort first,
% then give October and June 50% weight each. Do not let the unequal
% genotype counts in each cohort distort the visually displayed contrast.
[gc,kc]=findgroups(Mouse(:,{'genotype','cohort','context','phase_hz','amp_hz'}));
Cohort=kc;
Cohort.pac_z=splitapply(@mean,Mouse.pac_z,gc);
Cohort.n_mice=splitapply(@numel,Mouse.pac_z,gc);
writetable(Cohort,fullfile(output,'cohort_comodulograms_long.csv'));
[gg,kg]=findgroups(Cohort(:,{'genotype','context','phase_hz','amp_hz'}));
Group=kg;
Group.pac_z=splitapply(@mean,Cohort.pac_z,gg);
Group.n_cohorts=splitapply(@numel,Cohort.pac_z,gg);
assert(all(Group.n_cohorts==2),'Missing one of the two cohorts in group maps');
writetable(Group,fullfile(output,'group_comodulograms_long.csv'));

% Mouse-weighted ROI check against the published v2.0 baseline summaries.
[g,kg]=findgroups(Mouse(:,{'session_id','mouse_id','genotype','cohort','context'}));
Roi=kg;Roi.recomputed_roi_mean_z=splitapply(@mean,Mouse.pac_z,g);
RefMouse=readtable(refMouseFile,'TextType','string','VariableNamingRule','preserve');
check=innerjoin(Roi,RefMouse(:,{'session_id','context','morlet_roi_mean_z'}), ...
    'Keys',{'session_id','context'});
check.difference_z=check.recomputed_roi_mean_z-check.morlet_roi_mean_z;
writetable(check,fullfile(output,'roi_reference_comparison.csv'));
if mode=="locked"
    assert(height(check)==48,'ROI comparison missing mouse/context rows');
    assert(max(abs(check.difference_z))<0.001, ...
        'Recomputed ROI differs from frozen reference; see comparison CSV');
end

% Actual EEG-derived figures (maps are NOT pixel-wise significance tests).
fig_nrem(Mouse,figDir,C);
fig_rem(Mouse,figDir,C);
fig_mouse_rois(Roi,figDir);
fprintf('DONE. Baseline maps for %d mice. Files: %s\n', ...
    numel(unique(string(Mouse.mouse_id))),output);
fprintf('IMPORTANT: Group maps descriptive; inferential statistics use mouse-level ROI.\n');

%% ===== Local functions (same transform/normalization as run_20) =====
function c=cycles_from_filter(freqs,fs,halfWidth)
 c=zeros(size(freqs));
 for k=1:numel(freqs)
    f=freqs(k);lo=max(0.4,f-halfWidth);hi=min(fs/2-1,f+halfWidth);
    [b,a]=butter(2,[lo hi]/(fs/2),'bandpass');
    [hh,ff]=freqz(b,a,16384,fs);
    mag=abs(hh).^2;
    [peak,pidx]=max(mag);
    ixL=pidx;ixR=pidx;
    while ixL>1 && mag(ixL-1)>=peak/2,ixL=ixL-1;end
    while ixR<numel(ff) && mag(ixR+1)>=peak/2,ixR=ixR+1;end
    bw=ff(ixR)-ff(ixL);
    assert(isfinite(bw)&&bw>0,'Invalid FWHM at %g Hz',f);
    c(k)=min(25,max(2.5,2*sqrt(2*log(2))*f/bw));
 end
end
function [phase,amp]=morlet_components(x,fs,pf,af,ctx,C)
 x=double(x(:)).';x=x-mean(x);
 pad=min(numel(x)-1,round(C.padding_s*fs));
 xx=[fliplr(x(1:pad)) x fliplr(x(end-pad+1:end))];
 phase=zeros(numel(pf),numel(x));amp=zeros(numel(af),numel(x));
 if ctx=="REM",pc=C.cycles_pr;ac=C.cycles_ar;
 else,pc=C.cycles_pn;ac=C.cycles_an;end
 for j=1:numel(pf)
    z=morlet_conv(xx,pf(j),pc(j),fs,C.wavelet_support_sigma);
    phase(j,:)=angle(z(pad+(1:numel(x))));
 end
 for j=1:numel(af)
    z=morlet_conv(xx,af(j),ac(j),fs,C.wavelet_support_sigma);
    amp(j,:)=abs(z(pad+(1:numel(x))));
 end
end
function z=morlet_conv(x,f,cycles,fs,support)
 sigma=cycles/(2*pi*f);
 hw=ceil(support*sigma*fs);
 tw=(-hw:hw)/fs;
 w=exp(1i*2*pi*f*tw).*exp(-tw.^2/(2*sigma^2));
 w=w-mean(w);
 w=w/sqrt(sum(abs(w).^2));
 n=numel(x);nfft=2^nextpow2(n+numel(w)-1);
 ft=ifft(fft(x,nfft).*fft(w,nfft));
 z=ft(hw+(1:n));
end
function z=surrogate_zmap(phase,amp,shifts,C)
 n=size(phase,2);
 assert(n==size(amp,2)&&n>2*C.fs*C.min_shift_s,'Bad PAC support');
 Q=exp(1i*phase);
 raw=abs((amp*Q.')/n);
 sur=zeros([size(raw),numel(shifts)]);
 for k=1:numel(shifts)
    ap=circshift(amp,[0 shifts(k)]);
    sur(:,:,k)=abs((ap*Q.')/n);
 end
 mu=mean(sur,3);sd=std(sur,0,3);
 z=(raw-mu)./max(sd,1e-9);
end
function seed=window_seed(sid,context,bidx,ws,base)
 str=sprintf('%s|%s|%.0f|%.6f',char(sid),char(context),bidx,ws);
 h=double(base);
 for k=1:numel(str),h=mod(h*131+double(str(k)),2147483647);end
 seed=max(1,round(h));
end
function s=reclassify_short_wake(s,maxSeconds)
 s=s(:);d=diff([false;s==0;false]);a=find(d==1);b=find(d==-1)-1;
 for j=1:numel(a)
    if b(j)-a(j)+1<=maxSeconds,s(a(j):b(j))=15;end
 end
end
function T=state_runs(st)
 v=st(:);edges=[1;find(diff(v)~=0)+1;numel(v)+1];
 a=edges(1:end-1);b=edges(2:end)-1;
 T=table(a,b,v(a),'VariableNames',{'start_idx','end_idx','state'});
end
function [M,pf,af]=mean_map(Mouse,genotype,ctx)
 subset=Mouse(string(Mouse.genotype)==genotype & string(Mouse.context)==ctx,:);
 assert(~isempty(subset),'No observations for %s %s',char(genotype),char(ctx));
 pf=sort(unique(double(subset.phase_hz)));
 af=sort(unique(double(subset.amp_hz)));
 M=nan(numel(af),numel(pf));
 for i=1:numel(pf)
    for j=1:numel(af)
        q=double(subset.phase_hz)==pf(i) & double(subset.amp_hz)==af(j);
        v=double(subset.pac_z(q));
        assert(numel(v)==8,'Expected 8 mice for each genotype frequency cell');
        vt=double(subset.pac_z(q & string(subset.cohort)=="thesis"));
        vn=double(subset.pac_z(q & string(subset.cohort)=="new"));
        assert(~isempty(vt)&&~isempty(vn),'Missing cohort data');
        M(j,i)=0.5*(mean(vt)+mean(vn));
    end
 end
end
function M=contrast_map(Mouse,genotype,pf,af)
 ctx1=Mouse(string(Mouse.genotype)==genotype & string(Mouse.context)=="NREM_to_REM",:);
 ctx2=Mouse(string(Mouse.genotype)==genotype & string(Mouse.context)=="NREM_to_WAKE",:);
 M=zeros(numel(af),numel(pf));
 for cohort=["thesis","new"]
    mice=unique(string(ctx1.mouse_id(string(ctx1.cohort)==cohort)));
    assert(~isempty(mice),'Missing cohort in transition contrast');
    tmp=zeros(numel(af),numel(pf));
    for k=1:numel(mice)
        a=ctx1(string(ctx1.mouse_id)==mice(k),:);
        b=ctx2(string(ctx2.mouse_id)==mice(k),:);
        assert(~isempty(a)&&~isempty(b),'Missing transition context for %s',char(mice(k)));
        for i=1:numel(pf)
            for j=1:numel(af)
                qa=double(a.phase_hz)==pf(i) & double(a.amp_hz)==af(j);
                qb=double(b.phase_hz)==pf(i) & double(b.amp_hz)==af(j);
                tmp(j,i)=tmp(j,i)+double(a.pac_z(qa)-b.pac_z(qb))/numel(mice);
            end
        end
    end
    M=M+0.5*tmp;
 end
end
function lim=plot_limit(mats,minLim)
 a=[];
 for k=1:numel(mats),a=[a;abs(mats{k}(:))];end %#ok<AGROW>
 a=sort(a(isfinite(a)));
 if isempty(a),lim=minLim;return;end
 lim=max(minLim,ceil(10*a(max(1,ceil(0.98*numel(a)))))/10);
end
function draw_map(ax,pf,af,Z,ttl,lim)
 imagesc(ax,pf,af,Z);axis(ax,'xy');axis(ax,'tight');
 ax.XTick=pf(1:max(1,round(numel(pf)/5)):end);
 ax.YTick=af(1:max(1,round(numel(af)/6)):end);
 ax.FontName='Arial';ax.FontSize=10;ax.Box='off';
 clim(ax,[-lim lim]);colormap(ax,diverge(256));
 xlabel(ax,'Phase (Hz)');ylabel(ax,'Amplitude (Hz)');title(ax,ttl,'FontWeight','normal');
 cb=colorbar(ax);cb.Label.String='Mean PAC Z (or difference)';
end
function cmap=diverge(n)
 % Fixed blue-white-red divergence (neutral white = 0).
 s=linspace(0,1,n)';
 x=[0 0.5 1];
 cmap=[interp1(x,[0.15 0.97 0.73],s), ...
       interp1(x,[0.32 0.97 0.14],s), ...
       interp1(x,[0.67 0.97 0.16],s)];
end
function save_figure(f,d,name)
 exportgraphics(f,fullfile(d,[char(name) '.png']),'Resolution',300);
 exportgraphics(f,fullfile(d,[char(name) '.pdf']),'ContentType','vector');
 savefig(f,fullfile(d,[char(name) '.fig']));close(f);
end
function fig_nrem(Mouse,d,C)
 [a,pf,af]=mean_map(Mouse,"WT","NREM_to_REM");
 b=mean_map(Mouse,"WT","NREM_to_WAKE");
 c=mean_map(Mouse,"APP","NREM_to_REM");
 e=mean_map(Mouse,"APP","NREM_to_WAKE");
 wt=contrast_map(Mouse,"WT",pf,af);
 app=contrast_map(Mouse,"APP",pf,af);
 lim1=plot_limit({a,b,c,e},0.25);lim2=plot_limit({wt,app},0.15);
 f=figure('Color','w','Visible','off','Position',[40 40 1360 790]);
 tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
 for row=1:2
     if row==1,first=a;second=b;third=wt;gen='WT';
     else,first=c;second=e;third=app;gen='APP/PS1';end
     draw_map(nexttile(tl,(row-1)*3+1),pf,af,first,[gen ' | pre-REM'],lim1);
     draw_map(nexttile(tl,(row-1)*3+2),pf,af,second,[gen ' | pre-Wake'],lim1);
     draw_map(nexttile(tl,(row-1)*3+3),pf,af,third,[gen ' | pre-REM - pre-Wake'],lim2);
 end
 title(tl,'Baseline NREM slow-sigma PAC | cohort-balanced real Morlet maps');
 save_figure(f,d,'nrem_group_comodulograms');
 % Direct APP-WT maps and genotype difference-in-contrasts.
 f=figure('Color','w','Visible','off','Position',[40 40 1330 410]);
 tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
 mats={c-a,e-b,app-wt};ttl={'APP - WT | pre-REM','APP - WT | pre-Wake', ...
     'APP - WT | (pre-REM - pre-Wake)'};
 lim=plot_limit(mats,0.15);
 for j=1:3,draw_map(nexttile(tl),pf,af,mats{j},ttl{j},lim);end
 title(tl,'Descriptive cohort-balanced genotype differences (not pixelwise tests)');
 save_figure(f,d,'nrem_genotype_difference_comodulograms');
end
function fig_rem(Mouse,d,C) %#ok<INUSD>
 [wt,pf,af]=mean_map(Mouse,"WT","REM");app=mean_map(Mouse,"APP","REM");
 lim=plot_limit({wt,app},0.25);difflim=plot_limit({app-wt},0.15);
 f=figure('Color','w','Visible','off','Position',[40 40 1380 440]);
 tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
 draw_map(nexttile(tl),pf,af,wt,'WT | REM theta-gamma',lim);
 draw_map(nexttile(tl),pf,af,app,'APP/PS1 | REM theta-gamma',lim);
 draw_map(nexttile(tl),pf,af,app-wt,'APP - WT | REM',difflim);
 title(tl,'Baseline REM theta-gamma PAC | cohort-balanced real Morlet maps');
 save_figure(f,d,'rem_theta_gamma_group_comodulograms');
end
function fig_mouse_rois(Roi,d)
 f=figure('Color','w','Visible','off','Position',[40 40 1360 440]);
 tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
 contexts=["NREM_to_REM","NREM_to_WAKE","REM"];
 for i=1:3
    ax=nexttile(tl);hold(ax,'on');
    for j=1:height(Roi)
        if string(Roi.context(j))~=contexts(i),continue;end
        isAPP=string(Roi.genotype(j))=="APP";
        isNew=string(Roi.cohort(j))=="new";
        x=1+double(isAPP)+(double(isNew)-0.5)*0.18;
        if isAPP,col=[0.77 0.35 0.25];else,col=[0.19 0.41 0.70];end
        if isNew,mrk='s';else,mrk='o';end
        scatter(ax,x,double(Roi.recomputed_roi_mean_z(j)),34,col,mrk,'filled');
    end
    xlim(ax,[0.5 2.5]);xticks(ax,[1 2]);xticklabels(ax,{'WT','APP/PS1'});
    ylabel(ax,'Bout-balanced mean PAC Z');
    title(ax,char(contexts(i)),'Interpreter','none');
    ax.FontName='Arial';box(ax,'off');grid(ax,'on');
 end
 title(tl,'Mouse-level ROI summaries | circles: October, squares: June');
 save_figure(f,d,'mouse_level_roi_pac');
end
