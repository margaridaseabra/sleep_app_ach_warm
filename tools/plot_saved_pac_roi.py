#!/usr/bin/env python3
"""Plot REAL baseline Morlet PAC mouse-level ROI values already saved by v2.0.

No EEG recomputation, no guessed comodulogram cells, and no significance tests.
Requires: numpy, matplotlib (pandas NOT required).
"""
from __future__ import annotations
import argparse
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=(Path.home() / 'sleep_app_ach_warm' /
       'outputs/eeg/pac/baseline_morlet_bandwidth_v20/paired_pac_mouse.csv'))
    parser.add_argument('--out', type=Path, default=(Path.home() / 'sleep_app_ach_warm' /
       'outputs/eeg/pac/baseline_morlet_bandwidth_v20/paper_figures'))
    a = parser.parse_args()
    assert a.input.is_file(), f'Missing {a.input}'
    with a.input.open(newline='') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==48, f'Expected 16 mice × 3 contexts = 48 rows, found {len(rows)}'
    a.out.mkdir(parents=True, exist_ok=True)
    data={ (r['mouse_id'],r['context']):r for r in rows}
    mice={r['mouse_id']:(r['genotype'],r['cohort']) for r in rows}
    assert len(mice)==16 and len(data)==48
    measures={m:{k:float(data[(m,k)]['morlet_roi_mean_z']) for k in ['NREM_to_REM','NREM_to_WAKE','REM']}
              for m in mice}
    vals={m:[measures[m]['NREM_to_REM'],measures[m]['NREM_to_WAKE'],
             measures[m]['NREM_to_REM']-measures[m]['NREM_to_WAKE'],measures[m]['REM']]
          for m in mice}
    panels=[('Pre-REM NREM',0),('Pre-Wake NREM',1),('Pre-REM minus pre-Wake',2),('REM theta–gamma',3)]
    color={'WT':'#37689f','APP':'#b75a48'}
    fig,axs=plt.subplots(2,2,figsize=(10.5,7.6),layout='constrained')
    for ax,(name,j) in zip(axs.ravel(),panels):
        for geno,x in [('WT',1),('APP',2)]:
            v=[vals[m][j] for m,(g,c) in mice.items() if g==geno]
            for m,(g,c) in sorted(mice.items()):
                if g!=geno: continue
                y=vals[m][j]
                xx=x+(-0.13 if c=='thesis' else 0.13)
                ax.scatter(xx,y, s=42, marker=('o' if c=='thesis' else 's'),
                           facecolor=color[geno],edgecolor='white',linewidth=0.6,zorder=3)
            ax.plot([x-0.22,x+0.22],[np.mean(v)]*2,color='black',lw=2.0,zorder=4)
        ax.axhline(0,ls=':',color='0.6',lw=0.8)
        ax.set(xlim=(0.5,2.5),xticks=[1,2],xticklabels=['WT','APP/PS1'],
               ylabel='Bout-balanced mean PAC Z',title=name)
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y',alpha=0.16)
    fig.suptitle('Baseline Morlet PAC | 16 mice | actual data',fontsize=15,weight='bold')
    fig.text(0.5,-0.013,'Circles: October cohort     Squares: June cohort     Bars: group means; each point = one mouse',
             ha='center',fontsize=9)
    for ext in ['png','pdf','svg']:
        fig.savefig(a.out / ('baseline_morlet_mouse_level_ROIs.'+ext), dpi=300,bbox_inches='tight')
    plt.close(fig)
    print('Saved mouse-level PAC figure under',a.out)
    for geno in ['WT','APP']:
        print(geno,':', [round(np.mean([vals[m][j] for m,(g,c) in mice.items() if g==geno]),4)
                             for j in range(4)])

if __name__=='__main__':
    main()
