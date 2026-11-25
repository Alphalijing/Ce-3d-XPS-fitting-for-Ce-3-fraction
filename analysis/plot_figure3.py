"""Background constructions and the diagnostic region, from saved data only."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
d=pd.read_csv(ROOT/'data/validation/background/La_2at.csv')
out=ROOT/'results/reproduced_figures';out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'mathtext.fontset':'stix','pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(6.3,6.3),layout='constrained')
v=d[(d.energy_eV>=870)&(d.energy_eV<=895)]
high=v[v.high_mask]
ax.axvspan(high.energy_eV.min(),high.energy_eV.max(),color='#F1DE70',alpha=.28,label='High-confidence region')
ax.fill_between(v.energy_eV,v.noncrossing_Shirley/1000,v.converged_anchor_Shirley/1000,color='#9986BB',alpha=.28,label='Shirley background variation')
ax.plot(v.energy_eV,v.raw/1000,color='#222222',lw=.8,label='Measured spectrum')
ax.plot(v.energy_eV,v.noncrossing_Shirley/1000,color='#0072B2',lw=1.5,label='Noncrossing Shirley background')
ax.plot(v.energy_eV,v.converged_anchor_Shirley/1000,color='#D55E00',lw=1.4,label='Converged-anchor Shirley background')
ax.plot(v.energy_eV,v.linear/1000,color='#228833',lw=1.2,ls='--',label='Linear')
ax.set(xlim=(895,870),xlabel='Binding energy (eV)',ylabel=r'Intensity ($10^3$ counts s$^{-1}$)')
ax.tick_params(direction='in',top=True,right=True)
ax.legend(loc='upper left',fontsize=8.5,frameon=False)
for ext in ('png','pdf','svg'):fig.savefig(out/f'Figure3_background.{ext}',dpi=400)
plt.close(fig)
print('Saved background comparison to results/reproduced_figures/.')
