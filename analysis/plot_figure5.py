"""Synthetic recovery with conditional 95% posterior intervals; no fitting."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1]
m=pd.read_csv(ROOT/'data/synthetic/LA_MCMC_summary.csv')
f=pd.read_csv(ROOT/'data/synthetic/free_LA_summary.csv')
d=m.merge(f,on='case_id',suffixes=('_MCMC','_free'),validate='one_to_one')
out=ROOT/'results/reproduced_figures';out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'mathtext.fontset':'stix','pdf.fonttype':42,'svg.fonttype':'none'})
colors={.5:'#0072B2',1.:'#D55E00',1.5:'#009E73'}
fig,axes=plt.subplots(1,2,figsize=(12.4,6.4),sharex=True,sharey=True,layout='constrained')
for ax,panel,letter in zip(axes,('a','b'),('a','b')):
    ax.plot([-3,65],[-3,65],color='#999999',ls='--',lw=.9,zorder=0)
    for i,(sigma,col) in enumerate(colors.items()):
        for j,noise in enumerate((3,10)):
            a=d[(d.panel_MCMC==panel)&np.isclose(d.sigma_true_eV_MCMC,sigma)&np.isclose(d.noise_percent_MCMC,noise)].sort_values('Ce3_percent_true_window')
            x=a.Ce3_percent_true_window.to_numpy();offset=(i-1)*.8+(j-.5)*.32
            y=a.Ce3_percent_median.to_numpy();lo=a.Ce3_percent_q025.to_numpy();hi=a.Ce3_percent_q975.to_numpy()
            mark='o' if noise==3 else '^'
            ax.errorbar(x+offset,y,yerr=np.vstack([y-lo,hi-y]),fmt=mark,color=col,ms=5,capsize=2,lw=.8,zorder=3)
            ax.scatter(x+offset,a.Ce3_percent,marker=mark,s=32,facecolors='none',edgecolors=col,lw=1,zorder=2)
    ax.text(.03,.965,f'({letter})',transform=ax.transAxes,va='top',fontweight='bold',fontsize=13)
    ax.set(xlim=(-4,65),ylim=(-4,75),xticks=[0,15,30,60],xlabel=r'Preset Ce$^{3+}$ fraction (%)')
    ax.tick_params(direction='in',top=True,right=True)
    handles=[Line2D([],[],ls='none',marker='o',color='black',ms=5,label='Constrained Bayesian (95% CrI)'),
             Line2D([],[],ls='none',marker='o',mfc='none',mec='black',ms=5,label='Independent-parameter fit')]
    handles += [Line2D([],[],ls='none',marker='o',color=c,ms=5,label=rf'$\sigma$ = {s:g} eV') for s,c in colors.items()]
    handles += [Line2D([],[],ls='none',marker='o',color='#666666',label='3% noise'),Line2D([],[],ls='none',marker='^',color='#666666',label='10% noise')]
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(.02,.91),fontsize=8.1,frameon=False)
axes[0].set_ylabel(r'Fitted Ce$^{3+}$ fraction (%)')
for ext in ('png','pdf','svg'):fig.savefig(out/f'Figure5_synthetic.{ext}',dpi=400)
plt.close(fig)
print('Saved synthetic comparison to results/reproduced_figures/.')
