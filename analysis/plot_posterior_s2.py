"""One La composition per row, with a common Ce3+ percentage axis."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'results/reproduced_figures';out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'mathtext.fontset':'stix','pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(6,1,sharex=True,figsize=(7,10),layout='constrained')
for ax,at in zip(axes,(0,1,2,4,10,20)):
    with np.load(ROOT/f'data/experimental/la_series/posterior/LA_La_{at}at.npz') as z:draws=100*z['Ce3_fraction']
    q=np.quantile(draws,[.025,.5,.975])
    ax.hist(draws,bins=np.linspace(8,23,151),density=True,color='#0072B2',alpha=.7)
    ax.axvspan(q[0],q[2],color='#D55E00',alpha=.12)
    ax.axvline(q[1],color='#D55E00',lw=1.2)
    ax.text(.02,.90,f'{at} at% La',transform=ax.transAxes,va='top')
    ax.set_ylabel('Density');ax.tick_params(direction='in')
axes[-1].set(xlim=(8,23),xlabel=r'Ce$^{3+}$ fraction (%)')
for ext in ('png','pdf','svg'):fig.savefig(out/f'Figure_S2_LA_posteriors.{ext}',dpi=350)
plt.close(fig)
print('Saved six posterior distributions to results/reproduced_figures/.')
