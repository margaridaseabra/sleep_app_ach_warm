#!/usr/bin/env python3
"""Figure 3 ambient warming and PAC supplements from existing CSV results.

No EEG processing, no new inferential tests. Unit of plotted replication: mouse.
All numerical data come from previously computed analysis outputs.
"""
from pathlib import Path
import argparse, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

WT='#2879A9'; APP='#D46F50'; GREY='#697480'
G={'WT':WT,'APP':APP}; C={'thesis':'o','new':'D'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})

def need(p):
    if not p.is_file(): raise FileNotFoundError(f'Missing required input: {p}')
    return p

def save(fig,p):
    p.parent.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf','svg']:
        if ext=='png': fig.savefig(p.with_suffix('.png'),dpi=300,facecolor='white',bbox_inches='tight')
        else: fig.savefig(p.with_suffix('.'+ext),bbox_inches='tight',facecolor='white')
    plt.close(fig)

def groups_csv(df,value_col,outfile,percent=False):
    cols=[]
    for cohort,geno,label in [('thesis','WT','WT_Oct'),('thesis','APP','APP_Oct'),('new','WT','WT_June'),('new','APP','APP_June')]:
        q=df[(df.cohort==cohort)&(df.genotype==geno)].sort_values('mouse_id')
        vals=pd.to_numeric(q[value_col],errors='coerce').dropna().tolist()
        cols.append(pd.Series(vals,name=label))
    pd.concat(cols,axis=1).to_csv(outfile,index=False)

def miceplot(ax,df,metric,title,ylabel,stats=None,unit='',subtitle=''):
    d=df[df.metric==metric].copy()
    if d.empty: raise ValueError('No metric '+metric)
    xs={'WT':0,'APP':1}
    for genotype in ['WT','APP']:
        vals=[]
        for coh in ['thesis','new']:
            q=d[(d.genotype==genotype)&(d.cohort==coh)]
            n=len(q)
            if n==0: continue
            off={'thesis':-.075,'new':.075}[coh]
            x=np.full(n,xs[genotype]+off)
            y=q.difference_in_differences.to_numpy(float)
            ax.scatter(x,y,s=47,marker=C[coh],facecolor=G[genotype],edgecolor='white',linewidth=.6,alpha=.93,zorder=3)
            vals.extend(y[np.isfinite(y)])
        if vals:
            a=np.array(vals)
            mu=a.mean();sem=a.std(ddof=1)/np.sqrt(len(a)) if len(a)>1 else np.nan
            ax.errorbar([xs[genotype]],[mu],yerr=[sem],fmt='_',color='#20262d',capsize=6,elinewidth=1.8,markersize=18,zorder=5)
    ax.axhline(0,ls='--',color='#999da1',lw=1)
    ax.set_xticks([0,1],['WT','APP/PS1']);ax.set_xlim(-.45,1.45)
    ax.set_ylabel(ylabel);ax.set_title(title,loc='left',fontweight='bold',fontsize=11)
    if stats is not None:
        ax.text(.99,.99,f"n={int(stats.n_WT)} WT / {int(stats.n_APP)} APP\np={float(stats.p_welch):.3f} (uncorrected)",transform=ax.transAxes,fontsize=8,color='#555',ha='right',va='top',bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':1.2})

def rem_hr(ax,s):
    q=s[(s.time_bin=='ZT6_9')&(s.state=='REM') & (s.cohort_view.isin(['thesis','new']))].copy()
    if len(q)!=2: raise ValueError('Expected 2 cohort REM HR rows for ZT6_9')
    for i,cohort in enumerate(['thesis','new']):
        r=q[q.cohort_view==cohort].iloc[0]
        for gen,j in [('WT',0),('APP',1)]:
            stem=gen.lower(); p='wt' if gen=='WT' else 'app'
            y=i*3+j+.2
            center=r['HR_warming_vs_baseline_'+gen]
            low=r['CI_'+p+'_low']; high=r['CI_'+p+'_high']
            ax.errorbar(center,y,xerr=[[center-low],[high-center]],fmt='o',color=G[gen],ms=7,capsize=3)
        ax.text(0.98,i*3+1.52,f"{cohort}: interaction p={r.p_interaction:.3g}",transform=ax.get_yaxis_transform(),ha='right',fontsize=8,color=GREY)
    ax.axvline(1,color='#777',ls='--',lw=1)
    ax.set_xscale('log');ax.set_xlabel('Warming / baseline REM hazard ratio (95% CI)')
    ax.set_yticks([.2,1.2,3.2,4.2],['Oct WT','Oct APP','June WT','June APP'])
    ax.set_ylim(-.5,5.2);ax.set_xlim(.18,3.5)
    ax.set_title('F  REM bout termination hazard',loc='left',fontweight='bold',fontsize=11)
    ax.text(0.01,-.25,'Hazard ratio below 1 indicates longer REM bouts under warming.',transform=ax.transAxes,fontsize=8,color=GREY)

def timeline(ax):
    ax.set_title('A  Temperature manipulation protocol',loc='left',fontweight='bold',fontsize=11)
    ax.set_xlim(0,9);ax.set_ylim(0,1)
    for a,b,color,label in [(0,3,'#DCE7EE','Pre-warming'),(3,6,'#E9B887','Warming'),(6,9,'#E4E6E8','Washout')]:
        ax.add_patch(FancyBboxPatch((a,.35),b-a,.25,boxstyle='round,pad=0.015,rounding_size=.055',facecolor=color,edgecolor='#A3A8AD',linewidth=.6))
        ax.text((a+b)/2,.475,label,ha='center',va='center',fontsize=9)
    ax.set_xticks([0,3,6,9],['0','3','6','9'])
    ax.set_xlabel('Illustrative hours from recording start');ax.set_yticks([])
    for x in ['left','bottom','right','top']:ax.spines[x].set_visible(False)
    ax.text(0,.11,'Actual warming times vary by animal.',fontsize=8,color=GREY)

def km(df):
    d=df.sort_values('time_s')
    t=np.asarray(d.time_s,float);e=np.asarray(d.event,float)
    unique=np.unique(t[np.isfinite(t)])
    xx=[0.];yy=[1.];n=len(t);s=1.
    for tt in unique:
        at=(t==tt);de=np.sum(at & (e==1));at_risk=np.sum(t>=tt)
        if de and at_risk:
            s*=(1-de/at_risk)
            xx.append(float(tt));yy.append(float(s))
    return xx,yy

def km_figure(surv,p):
    fig,axarr=plt.subplots(2,2,figsize=(12.6,9.3),sharey=True)
    for i,coh in enumerate(['thesis','new']):
      for j,state in enumerate(['NREM','REM']):
        ax=axarr[i,j]
        for gt in ['WT','APP']:
          for condition in ['baseline','warming']:
            df=surv[(surv.cohort==coh)&(surv.genotype==gt)&(surv.time_bin=='ZT6_9')&(surv.state==state)&(surv.condition==condition)]
            if len(df):
                x,y=km(df)
                ax.step(x,y,where='post',color=G[gt],ls='-' if condition=='warming' else '--',lw=1.65,label=f'{gt}, {condition}')
        ax.set_title(f'{"October" if coh=="thesis" else "June"} · {state}',fontweight='bold')
        ax.set_xlabel('Bout duration (s)');ax.set_ylabel('Fraction of bouts surviving')
        ax.set_xlim(0,220 if state=='REM' else 320);ax.set_ylim(0,1);ax.grid(alpha=.14)
    axarr[0,0].legend(ncol=2,fontsize=8,loc='upper right')
    fig.suptitle('Descriptive bout survival, ZT6–9 · cohorts kept separate',fontsize=14,fontweight='bold',y=.995)
    fig.text(.02,.014,'Bouts are not independent mice. Statistical interpretation requires mouse-clustered models; see original Cox tables.',fontsize=8,color=GREY)
    fig.tight_layout(rect=[0,.03,1,.97]);save(fig,p)

def pac_supp(root,outf,outpr):
    p= root/'outputs/eeg/pac/baseline_morlet_bandwidth_v20/paired_pac_mouse.csv'
    pac=pd.read_csv(need(p))
    for context,short in [('NREM_to_REM','NREM_preREM'),('NREM_to_WAKE','NREM_preWake'),('REM','REM_theta_gamma')]:
        d=pac[pac.context==context]
        groups_csv(d,'morlet_roi_mean_z',outpr/f'Prism_baseline_{short}_Morlet_Z.csv')
    wide=pac.pivot_table(index=['mouse_id','cohort','genotype'],columns='context',values='morlet_roi_mean_z').reset_index()
    wide['preREM_minus_preWake']=wide.NREM_to_REM-wide.NREM_to_WAKE
    wide.to_csv(outpr/'PAC_baseline_mouse_level_long.csv',index=False)
    groups_csv(wide,'preREM_minus_preWake',outpr/'Prism_baseline_NREM_transition_contrast_Morlet_Z.csv')
    fig,axs=plt.subplots(1,3,figsize=(12.2,4.3))
    for ax,(column,title) in zip(axs,[('NREM_to_REM','A  Pre-REM NREM slow–sigma'),('preREM_minus_preWake','B  NREM transition contrast'),('REM','C  REM theta–gamma')]):
        d=wide[['mouse_id','cohort','genotype',column]].rename(columns={column:'difference_in_differences'}).copy()
        d['metric']='roi'
        miceplot(ax,d,'roi',title,'Morlet PAC Z')
    axs[0].text(0,-.28,'Phase 1–4 Hz / amplitude 9–16 Hz',transform=axs[0].transAxes,fontsize=8,color=GREY)
    axs[2].text(0,-.28,'Phase 4–10 Hz / amplitude 30–80 Hz',transform=axs[2].transAxes,fontsize=8,color=GREY)
    fig.text(.015,.01,'Mouse-level descriptive values; no pixelwise inference. Circle: October, diamond: June.',fontsize=9,color=GREY)
    fig.tight_layout(rect=[0,.05,1,1]);save(fig,outf/'Fig1S_PAC_mouse_level')
    real=root/'outputs/eeg/pac/real_comodulograms_v22_fine/figures'
    for name in ['nrem_group_comodulograms','rem_theta_gamma_group_comodulograms']:
        for ext in ['png','pdf']:
            src=real/(name+'.'+ext)
            if src.is_file(): shutil.copy2(src,outf/(f'Fig1S_real_{name}.{ext}'))

def build(args):
    root=Path(args.repo).expanduser().resolve()
    base=root/'outputs/sleep_architecture/warming_complete'
    d=pd.read_csv(need(base/'warming_difference_in_differences.csv'))
    stats=pd.read_csv(need(base/'warming_DiD_genotype_stats.csv'))
    survstat=pd.read_csv(need(base/'survival/warming_baseline_vs_warming_cox_interaction.csv'))
    surv=pd.read_csv(need(base/'survival/warming_survival_dataset.csv'))
    outfig=root/'results/figures/paper/figure3_warming';outpr=root/'results/prism/paper/figure3_warming'
    pacfig=root/'results/figures/paper/figure1_baseline/pac_extended';pacpr=root/'results/prism/paper/figure1_baseline/pac_extended'
    for p in [outfig,outpr,pacfig,pacpr]:p.mkdir(parents=True,exist_ok=True)
    warm=d[d.time_bin=='ZT6_9'].copy()
    if len(warm[warm.metric=='percent_time'].mouse_id.unique())!=16:
        raise AssertionError('Expected 16 matched mice in primary ZT6–9 analysis')
    warm.to_csv(outpr/'warming_mouse_level_DiD_all_metrics_LONG.csv',index=False)
    stats.to_csv(outpr/'original_warming_genotype_stats_reference.csv',index=False)
    survstat.to_csv(outpr/'original_warming_survival_interaction_reference.csv',index=False)
    configs=[('state_percent','REM','percent_time','REM_percent_DiD_pp'),('bout_metric','REM','mean_bout_duration_s','REM_duration_DiD_s'),('state_percent','NREM','percent_time','NREM_percent_DiD_pp'),('state_percent','MA','percent_time','MA_percent_DiD_pp')]
    fig,axs=plt.subplots(2,3,figsize=(15.5,9.3))
    timeline(axs[0,0])
    titles=['B  REM time','C  REM bout duration','D  NREM time','E  Microarousal time']
    for a,(dom,st,metric,name),title in zip([axs[0,1],axs[0,2],axs[1,0],axs[1,1]],configs,titles):
        sub=warm[(warm.domain==dom)&(warm.state==st)&(warm.metric==metric)]
        if len(sub)!=16: raise AssertionError(f'Expected 16 values for {name}, got {len(sub)}')
        sub[['mouse_id','genotype','cohort','baseline_pre','baseline_target','warming_pre','warming_target','difference_in_differences']].to_csv(outpr/f'{name}_mouse_level_long.csv',index=False)
        groups_csv(sub,'difference_in_differences',outpr/f'Prism_warming_{name}.csv')
        s=stats[(stats.cohort_view=='all')&(stats.time_bin=='ZT6_9')&(stats.domain==dom)&(stats.state==st)&(stats.metric==metric)]
        miceplot(a,sub,metric,title,'Matched baseline-controlled change\n'+('percentage points' if name.endswith('pp') else 'seconds'),s.iloc[0] if len(s) else None)
    rem_hr(axs[1,2],survstat)
    fig.suptitle('Figure 3 · Warming responses differ by genotype, cohort and sleep endpoint',fontsize=15,fontweight='bold')
    fig.text(.01,.008,'Mouse-level difference-in-differences: (warming active − pre) − (matched baseline active − pre). Error bars are SEM. Circle: October, diamond: June. ',fontsize=8,color=GREY)
    fig.tight_layout(rect=[0,.026,1,.965],h_pad=2.5,w_pad=2.2)
    save(fig,outfig/'Figure3_warming_data_six_panel_draft')
    km_figure(surv,outfig/'Fig3S_NREM_REM_survival_cohorts')
    # PAC supplement for warming (and per-mouse Prism-ready data)
    pacpath=root/'outputs/eeg/pac/intervention_hilbert_v15/matched_review/paired_intervention_mouse.csv'
    if pacpath.is_file():
        pac=pd.read_csv(pacpath)
        w=pac[(pac.intervention=='warming')&(pac.n_pre>=3)&(pac.n_active>=3)&(pac.n_matched_baseline_pre>=3)&(pac.n_matched_baseline_active>=3)].copy()
        w.to_csv(outpr/'warming_Hilbert_PAC_matched_mouse_LONG.csv',index=False)
        f,a=plt.subplots(1,3,figsize=(12.5,4.1))
        for ax,(context,title) in zip(a,[('NREM_to_REM','Pre-REM slow–sigma'),('NREM_to_WAKE','Pre-Wake slow–sigma'),('REM','REM theta–gamma')]):
            q=w[w.context==context].copy(); q=q.rename(columns={'roi_mean_z_DiD':'difference_in_differences'});q['metric']='roi'
            miceplot(ax,q,'roi',title,'Matched PAC response (Z)')
            groups_csv(q,'difference_in_differences',outpr/f'Prism_warming_Hilbert_PAC_{context}_DiD_Z.csv')
        f.suptitle('Exploratory warming PAC (Hilbert; no confirmed genotype-specific rescue)',fontweight='bold')
        f.tight_layout(rect=[0,.01,1,.95]);save(f,outfig/'Fig3S_warming_PAC_mouse_level')
    pac_supp(root,pacfig,pacpr)
    print('Figure 3 plots:',outfig);print('Figure 3 Prism:',outpr)
    print('Figure 1 PAC plots:',pacfig);print('Figure 1 PAC Prism:',pacpr)
    print('Original baseline, warming, spindle, and PAC outputs unchanged.')

if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--repo',default='.');build(pa.parse_args())
