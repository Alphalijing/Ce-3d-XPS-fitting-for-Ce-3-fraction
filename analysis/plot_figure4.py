"""Rearrange the archived Figure 4 data; no fitting or new estimates."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/reproduced_figures'
OUT.mkdir(parents=True,exist_ok=True)
TABLE=ROOT/'data/experimental/la_series/Figure4_composition_panels_data.csv'
LOCAL=ROOT/'data/experimental/la_series/Figure4_local_fit_data.csv'
COLORS={'Voigt':'#0072B2','LA':'#D55E00','DS':'#009E73'}
LABELS={'Voigt':'Voigt','LA':'LA-type','DS':'DS'}
MARKERS={'Voigt':'o','LA':'s','DS':'^'}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
before={str(p):sha(p) for p in (TABLE,LOCAL)}
d=pd.read_csv(TABLE)
curves=pd.read_csv(LOCAL)
xs=d.La_at_percent.to_numpy()
assert np.array_equal(xs,np.array([0,1,2,4,10,20]))
energy=curves.energy_eV.to_numpy()
raw=curves.raw.to_numpy()
background=curves.Voigt_fixed_background.to_numpy()
for model in COLORS:
    assert np.allclose(curves[f'{model}_fixed_background'],background,rtol=0,atol=1e-8)

plt.rcParams.update({'font.family':'Times New Roman','font.size':11,
    'mathtext.fontset':'stix','axes.linewidth':.9,'axes.labelsize':12,
    'xtick.direction':'in','ytick.direction':'in','xtick.labelsize':10,
    'ytick.labelsize':10,'legend.fontsize':11,'legend.frameon':False,
    'savefig.facecolor':'white','svg.fonttype':'none','pdf.fonttype':42})

fig=plt.figure(figsize=(11,7.6))
# Both blocks have a physical height of 6.156 inches. The left block's
# width is 3.08 inches, giving a 1.999:1 height/width ratio.
left=fig.add_gridspec(2,1,left=.10,right=.38,bottom=.12,top=.93,
    height_ratios=[3.8,1],hspace=.07)
fit=fig.add_subplot(left[0])
residual=fig.add_subplot(left[1],sharex=fit)
right=fig.add_gridspec(3,1,left=.50,right=.89,bottom=.12,top=.93,hspace=.065)
b=fig.add_subplot(right[0])
c=fig.add_subplot(right[1],sharex=b)
g=fig.add_subplot(right[2],sharex=b)

# (a) retains the measured leading edge, all fitted profiles and residuals.
fit.scatter(energy,(raw-background)/1000,s=9,color='#222222',linewidths=0,zorder=5)
net=[(raw-background)/1000]
local_RMSE={}
for model in COLORS:
    total=curves[f'{model}_fixed_total'].to_numpy()
    fit.plot(energy,(total-background)/1000,color=COLORS[model],lw=1.5,zorder=4 if model=='LA' else 3)
    residual.plot(energy,(raw-total)/1000,color=COLORS[model],lw=1.05)
    net.append((total-background)/1000)
    rmse=float(np.sqrt(np.mean((raw-total)**2)))
    assert np.isclose(rmse,d.loc[d.La_at_percent.eq(2),f'{model}_high_RMSE'].iloc[0],rtol=1e-12)
    local_RMSE[model]=rmse
net=np.concatenate(net)
fit.set(xlim=(float(energy.max()),float(energy.min())),
    ylim=(min(-.8,float(net.min())-.35),float(net.max())+1.5),
    ylabel=r'Net intensity ($10^3$ counts s$^{-1}$)')
fit.tick_params(axis='x',labelbottom=False)
fit.text(.045,.965,'(a)',transform=fit.transAxes,va='top',fontweight='bold',fontsize=13)
fit.text(.045,.91,'2 at% La',transform=fit.transAxes,va='top',fontsize=11)
residual.axhline(0,color='#888888',ls=':',lw=.9,zorder=0)
residual.set(xlabel='Binding energy (eV)',ylabel='Residual',
    ylim=(-3.1,2.1),xticks=[878,876,874,872,870],yticks=[-2,0,2])
residual.yaxis.label.set_size(11)

# (b,c,d) share a single, numerical La-content axis.
for model in COLORS:
    style=dict(color=COLORS[model],marker=MARKERS[model],markersize=5.1,
        linewidth=1.3,markeredgewidth=.9)
    b.plot(xs,d[f'{model}_high_RMSE'],**style)
    c.plot(xs,d[f'{model}_Ce3_window_percent'],**style)
    g.plot(xs,d[f'{model}_sigma_eV'],**style)

# Equal-tail 95% credible intervals for MCMC quantities only; no DS intervals.
for model in ('Voigt','LA'):
    for ax,key,low,high in ((c,'Ce3_window_percent','Ce3_low','Ce3_high'),
                            (g,'sigma_eV','sigma_low','sigma_high')):
        values=d[f'{model}_{key}'].to_numpy()
        errors=np.vstack([values-d[f'{model}_{low}'].to_numpy(),d[f'{model}_{high}'].to_numpy()-values])
        ax.errorbar(xs,values,yerr=errors,fmt='none',ecolor=COLORS[model],
            elinewidth=.8,capsize=2,capthick=.8,alpha=.62,zorder=1)
all_rmse=d[[m+'_high_RMSE' for m in COLORS]].to_numpy()
b_low=min(450,float(50*np.floor(np.min(all_rmse)/50)))
b_high=max(1250,float(50*np.ceil(np.max(all_rmse)/50)+50))
b.set(ylim=(b_low,b_high),yticks=np.arange(500,b_high+1,250),
    ylabel='High-confidence RMSE\n'+r'(counts s$^{-1}$)')
b.yaxis.set_label_position('right');b.yaxis.tick_right()
b.yaxis.labelpad=10
c.set(ylim=(min(10,float(np.floor(d[['Voigt_Ce3_low','LA_Ce3_low']].min().min())-1)),24),yticks=[10,12,14,16,18,20,22,24],ylabel=r'Ce$^{3+}$ fraction (%)')
c.axvline(10,color='#929292',ls=':',lw=.9,zorder=0)
catalysis=c.twinx()
catalysis.plot(xs,d.DMC_productivity_mmol_g,color='#666666',marker='D',
    markersize=4.6,mfc='white',mec='#666666',linewidth=1.15,ls='--')
catalysis.set(ylim=(0,16),yticks=[0,4,8,12,16],
    ylabel=r'DMC productivity (mmol g$^{-1}$)')
catalysis.yaxis.label.set_color('#555555')
catalysis.tick_params(axis='y',labelcolor='#555555')
catalysis.yaxis.label.set_size(11)
catalysis.yaxis.labelpad=10
g.set(ylim=(min(.59,float(d[['Voigt_sigma_low','LA_sigma_low']].min().min())-.025),max(.94,float(d[['Voigt_sigma_high','LA_sigma_high']].max().max())+.025)),yticks=[.6,.7,.8,.9],ylabel=r'Gaussian $\sigma$ (eV)',
    xlabel='Nominal La content (at%)',xlim=(-.7,20.7),xticks=xs)
for ax in (b,c):ax.tick_params(axis='x',labelbottom=False)
for ax,letter in ((b,'b'),(c,'c'),(g,'d')):
    ax.text(.025,.965,f'({letter})',transform=ax.transAxes,va='top',
        fontweight='bold',fontsize=13)

# One shared legend for the entire figure.
handles=[Line2D([],[],color=COLORS[m],marker=MARKERS[m],lw=1.5,ms=5.2,
    label=LABELS[m]) for m in COLORS]
handles.extend([
    Line2D([],[],color='#222222',marker='o',ls='none',ms=3.4,label='Experiment'),
    Line2D([],[],color='#666666',marker='D',mfc='white',lw=1.15,ls='--',ms=4.6,label='DMC productivity')])
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.997),
    ncols=5,columnspacing=1.65,handletextpad=.55,handlelength=2.0,borderaxespad=.15)

stem=OUT/'Figure4_MCMC_Voigt_LA_DS_reference'
for ext in ('png','pdf','svg'):
    fig.savefig(stem.with_suffix('.'+ext),dpi=600,facecolor='white')
fig.savefig(OUT/'Figure4_MCMC_Voigt_LA_DS_reference_preview.png',dpi=180,facecolor='white')
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
# Record extents as an objective check for cropped labels.
issues=[]
visible_text=list(fig.texts)
for ax in fig.axes:
    visible_text.extend([ax.xaxis.label,ax.yaxis.label,*ax.texts])
    lo,hi=sorted(ax.get_xlim())
    visible_text.extend(t for v,t in zip(ax.get_xticks(),ax.get_xticklabels()) if lo<=v<=hi)
    lo,hi=sorted(ax.get_ylim())
    visible_text.extend(t for v,t in zip(ax.get_yticks(),ax.get_yticklabels()) if lo<=v<=hi)
for artist in visible_text:
    if not artist.get_visible() or not artist.get_text():continue
    bounds=artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    if bounds.x0<-.001 or bounds.y0<-.001 or bounds.x1>1.001 or bounds.y1>1.001:
        issues.append({'text':artist.get_text(),'figure_bounds':list(bounds.bounds)})
assert not issues,issues
plt.close(fig)
assert before=={str(p):sha(p) for p in (TABLE,LOCAL)}
manifest={'source_files_sha256':before,'fitting_performed':False,'data_modified':False,'estimators':{'LA':'MCMC posterior median','Voigt':'MCMC posterior median','DS':'retained least-squares point estimate'},'errorbars':'95% equal-tail credible intervals for LA and Voigt only',
    'figure_size_inches':[11,7.6],
    'panel_a_block_height_width_ratio':(.93-.12)*7.6/((.38-.10)*11),
    'left_and_right_block_height_inches':(.93-.12)*7.6,
    'shared_x_panels':['b','c','d'],'single_model_legend':True,
    'b_axis_position':'right','b_y_limits':[b_low,b_high],'DMC_y_limits':[0,16],
    'representative_sample_La_at_percent':2,'local_RMSE_CPS':local_RMSE,
    'cropped_text_checks':issues,
    'outputs':{p.name:sha(p) for p in OUT.glob('Figure4_MCMC_Voigt_LA_DS_reference*')}}
(OUT/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'figure':str(stem),'a_height_width_ratio':manifest['panel_a_block_height_width_ratio'],'data_unchanged':True},ensure_ascii=False))
