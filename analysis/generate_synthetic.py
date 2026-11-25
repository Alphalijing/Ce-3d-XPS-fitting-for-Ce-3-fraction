"""LA-theory synthetic spectra with a separately added cumulative loss tail.

This script generates data only. It never starts a fit or edits released code.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid

SOURCE=Path(__file__).resolve().parents[1]
import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--outdir',type=Path,default=SOURCE/'results/regenerated_synthetic')
_args=_parser.parse_args()
DATA=_args.outdir.resolve()
WORK=DATA
sys.path.insert(0,str(SOURCE))
from ce3d_xps_model import Ce3dModel, MAIN_WEIGHTS

GAMMAS=[.6436838494848398,.12294174414630221,.1572645154473184,.9943953028407514,.2681632412537586]
SEED=20261116


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_bytes(path,content):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists() and path.read_bytes()!=content:
        raise RuntimeError(f'Refusing to replace different existing synthetic data: {path}')
    path.write_bytes(content)


def save_json(path,obj):
    save_bytes(path,json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False).encode('utf-8'))


def save_csv(path,frame):
    save_bytes(path,frame.to_csv(index=False,float_format='%.17g').encode('utf-8-sig'))


def label(v):return f'{v:g}'.replace('.','p')


def generate():
    energy=np.round(np.arange(870.58,930.58+.001,.1),2)
    keep=energy<=925
    efit=energy[keep]
    full=Ce3dModel(energy,np.zeros_like(energy),lineshape='LA')
    fit=Ce3dModel(efit,np.zeros_like(efit),lineshape='LA')
    rows=[];checks=[]
    for batch in ('background_free','cumulative_scattering'):
        for x in (0.,.15,.3,.6):
            for sigma in (.5,1.,1.5):
                theta=np.array([0.,sigma,x,1.,*GAMMAS,0.,0.])
                profiles=full._unit_profiles(theta[None,:])[0]
                profiles/=np.trapz(profiles[keep],efit,axis=0)[None,:]
                amps=MAIN_WEIGHTS*np.where(np.arange(10)<4,x,1-x)
                unit_components=profiles*amps[None,:]
                unit=unit_components.sum(axis=1)
                np.testing.assert_allclose(unit[keep],fit.predict(theta)[0],rtol=2e-13,atol=1e-14)
                for noise in (.03,.10):
                    case=f'LA_Ce3_{round(100*x):03d}pct_sigma_{label(sigma)}eV_poisson_{round(100*noise):02d}pct'
                    case_seed=SEED+len(rows)
                    scale=(1/noise**2)/unit.max()
                    components=unit_components*scale;photo=components.sum(axis=1)
                    cumulative=cumulative_trapezoid(photo,energy,initial=0.)
                    tail_shape=cumulative/cumulative.max()*photo.max()
                    b0=0. if batch=='background_free' else .86*photo.max()
                    beta_tail=0. if batch=='background_free' else .84
                    tail=beta_tail*tail_shape;bg=b0+tail
                    mu=photo+bg
                    raw=np.random.default_rng(case_seed).poisson(mu).astype(float)
                    truth_areas=np.trapz(components[keep],efit,axis=0)
                    actual_x=truth_areas[:4].sum()/truth_areas.sum()
                    np.testing.assert_allclose(actual_x,x,rtol=0,atol=2e-15)
                    np.testing.assert_allclose(truth_areas.sum(),scale,rtol=2e-14)
                    assert tail[0]==0 and np.all(np.diff(tail)>=-1e-10)
                    out=DATA/batch
                    inp=out/f'{case}.csv'
                    clean=out/f'{case}_truth_curves.csv'
                    comp=out/f'{case}_components.csv'
                    save_csv(inp,pd.DataFrame({'Binding Energy (eV)':energy[::-1],'Intensity':raw[::-1],'True_background':bg[::-1]}))
                    save_csv(clean,pd.DataFrame({'energy_eV':energy,'photo':photo,'tail_shape_S':tail_shape,
                        'tail_beta_S':tail,'baseline_B0':np.full_like(energy,b0),'true_background':bg,
                        'raw_expectation':mu,'raw_noisy':raw,'in_fit_window':keep.astype(int)}))
                    cframe=pd.DataFrame({'energy_eV':energy})
                    for j,name in enumerate(full.component_names):cframe[name]=components[:,j]
                    save_csv(comp,cframe)
                    row=dict(batch=batch,case=case,file=str(inp),clean_reference_file=str(clean),components_file=str(comp),
                        Ce3_fraction_true=x,Ce3_percent_true=100*x,Ce3_percent_true_window=float(100*actual_x),
                        sigma_ce_true=sigma,poisson_relative_noise_at_photo_peak=noise,noise_percent=100*noise,
                        A_ce_total_true_window=float(scale),photo_peak_counts=float(photo.max()),
                        B0=float(b0),beta_tail=beta_tail,B0_to_photo_peak=0. if batch=='background_free' else .86,
                        LA_alpha=.9,LA_beta=2.,g3_1=GAMMAS[0],g3_2=GAMMAS[1],g4_1_LA_scale=GAMMAS[2],
                        g4_2=GAMMAS[3],g4_3=GAMMAS[4],seed=case_seed,
                        original_points=len(energy),fit_points=int(keep.sum()),energy_min_eV=float(energy[0]),
                        energy_max_eV=float(energy[-1]),fit_energy_max_eV=float(efit[-1]),
                        input_sha256=sha(inp),truth_curves_sha256=sha(clean),components_sha256=sha(comp))
                    save_json(out/f'{case}_truth.json',dict(**row,component_names=list(full.component_names),
                        component_areas_in_fit_window=truth_areas.tolist(),
                        raw_distribution='independent Poisson draws with mean raw_expectation',
                        background_formula='B0 + beta_tail * cumulative_trapezoid(photo)/max(cumulative_trapezoid(photo)) * max(photo)',
                        interval_normalization='each theoretical peak has unit area over 870.58–924.98 eV before assigning its component area'))
                    rows.append(row)
                    t=theta.copy();t[3]=scale
                    max_difference=float(np.max(abs(photo[keep]-fit.predict(t)[0])))
                    checks.append(dict(batch=batch,case=case,true_x_area_error=float(actual_x-x),
                        max_clean_difference_from_released_model=max_difference,
                        true_background_crosses_noisy_raw_points=int(np.sum(bg>raw)),
                        tail_at_Emin=float(tail[0]),tail_is_monotone=True,
                        raw_counts_nonnegative_integers=bool(np.all((raw>=0)&(raw==np.round(raw))))))
    save_csv(DATA/'synthetic_LA_ground_truth.csv',pd.DataFrame(rows))
    save_csv(WORK/'synthetic_LA_ground_truth.csv',pd.DataFrame(rows))
    meta=dict(description='Matched LA-theory synthetic benchmark; not independent validation of a microscopic lineshape.',
        spectra=48,batches={'background_free':24,'cumulative_scattering':24},
        theory_source=str(SOURCE),theory_version='2.0.1',
        source_sha256={n:sha(SOURCE/n) for n in ['ce3d_xps_model.py','ce3d_io.py','shirley_v7.py']},
        asymmetric_peaks=['Ce4+ v','Ce4+ u'],other_eight_peaks='Voigt',LA_alpha=.9,LA_beta=2.,
        g4_1_definition='Powered-Lorentz LA kernel scale; not Lorentzian HWHM',
        Ce3_fraction_definition='CeIII component areas / all component areas within 870.58–924.98 eV',
        full_energy_grid_eV=[870.58,930.58,.1],fit_energy_grid_eV=[870.58,924.98,.1],
        cumulative_integral='cumulative_trapezoid on ascending binding energy with initial=0',
        decay_kernel=None,beta_tail=.84,B0_to_photo_peak=.86,
        random_seeds=[SEED,SEED+47],noise='Poisson; nominal noise parameter scales photo peak to 1/noise^2 counts',
        note='Noise in the raw spectrum with background is not a uniform fixed relative percentage. Do not clip raw values to the true background.',
        fitting_performed=False,legacy_Voigt_truth_results_reused=False)
    save_json(DATA/'synthetic_LA_metadata.json',meta);save_json(WORK/'synthetic_LA_metadata.json',meta)
    save_json(WORK/'generation_checks.json',dict(cases=checks,
        max_area_fraction_error=max(abs(c['true_x_area_error']) for c in checks),
        max_clean_model_difference=max(c['max_clean_difference_from_released_model'] for c in checks),
        all_passed=True,sampling_performed=False))
    return pd.DataFrame(rows)


def preview(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'mathtext.fontset':'stix','svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(11,8.4),sharex=True,sharey='row',layout='constrained')
    for i,batch in enumerate(('background_free','cumulative_scattering')):
        selected=rows[(rows.batch==batch)&rows.Ce3_percent_true.isin([0.,60.])&(rows.sigma_ce_true==1.)&(rows.noise_percent==3.)]
        row_max=max(float(pd.read_csv(p).raw_noisy.max()) for p in selected.clean_reference_file)
        for j,x in enumerate((0.,60.)):
            r=rows[(rows.batch==batch)&(rows.Ce3_percent_true==x)&(rows.sigma_ce_true==1.)&(rows.noise_percent==3.)].iloc[0]
            d=pd.read_csv(r.clean_reference_file);ax=axes[i,j]
            ax.plot(d.energy_eV,d.raw_noisy,color='#444444',lw=.65,label='Simulated counts')
            ax.plot(d.energy_eV,d.raw_expectation,color='#0072B2',lw=1.1,label='Expected intensity')
            ax.plot(d.energy_eV,d.true_background,color='#D55E00',lw=1.2,ls='--',label='Generation background')
            ax.axvline(925,color='#BBBBBB',lw=.7,ls=':')
            ax.set_title(f'{"Background-free" if i==0 else "Cumulative scattering"} | preset Ce'+r'$^{3+}$'+f' = {x:g}%')
            ax.set_xlim(931,870);ax.set_ylim(-40,1.08*row_max)
            ax.text(.02,.96,f'({chr(97+i*2+j)})',transform=ax.transAxes,va='top',fontweight='bold')
            if i==1:ax.set_xlabel('Binding energy (eV)')
            if j==0:ax.set_ylabel('Simulated counts')
    axes[0,0].legend(loc='upper left',bbox_to_anchor=(.02,.88),fontsize=8,frameon=False)
    fig.suptitle(r'LA-theory synthetic spectra: $\sigma=1.0$ eV; nominal photo-peak noise 3%',fontsize=13)
    for ext in ('png','svg','pdf'):fig.savefig(WORK/f'LA_synthetic_generation_preview.{ext}',dpi=220)
    plt.close(fig)


if __name__=='__main__':
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    rows=generate();preview(rows)
    print(f'Generated {len(rows)} spectra in {DATA}',flush=True)
