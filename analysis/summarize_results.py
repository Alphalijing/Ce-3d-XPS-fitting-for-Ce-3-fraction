"""Recalculate publication statistics from distributed tables, without fitting."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]

def calculate():
    d=pd.read_csv(ROOT/'data/experimental/la_series/Figure4_composition_panels_data.csv')
    ext=pd.read_csv(ROOT/'data/experimental/additional_19/summary.csv')
    m=pd.read_csv(ROOT/'data/synthetic/LA_MCMC_summary.csv')
    f=pd.read_csv(ROOT/'data/synthetic/free_LA_summary.csv')
    paired=m.merge(f,on='case_id',suffixes=('_MCMC','_free'),validate='one_to_one')
    me=paired.Ce3_percent_median-paired.Ce3_percent_true_window
    fe=paired.Ce3_percent-paired.Ce3_percent_true_window
    lower=paired.RMSE_free<paired.RMSE_MCMC
    larger=fe.abs()>me.abs()
    records=[]
    for panel in ('a','b'):
        for noise in (3,10):
            mask=(paired.panel_MCMC==panel)&np.isclose(paired.noise_percent_MCMC,noise)
            for name,err in (('MCMC',me),('independent_LA_LS',fe)):
                values=err[mask].to_numpy()
                records.append(dict(panel=panel,noise_percent=noise,method=name,n=len(values),
                    MAE_pp=float(np.mean(abs(values))),mean_signed_error_pp=float(np.mean(values)),
                    min_error_pp=float(values.min()),max_error_pp=float(values.max())))
    joint=pd.read_csv(ROOT/'data/validation/joint_background/all_fixed_and_joint_fits.csv')
    j=joint[joint.background=='joint'].merge(joint[joint.background=='fixed'],on=['sample','shape'],suffixes=('_joint','_fixed'))
    ok=j.success_joint.astype(str).str.lower().eq('true')
    accepted=j[ok]
    boundary=accepted.bounds_joint.fillna('').str.contains('anchor_')
    gains=100*(1-d.LA_high_RMSE/d.Voigt_high_RMSE)
    return dict(
        main_LA_vs_Voigt_RMSE_reduction_percent=gains.tolist(),
        main_LA_better_than_DS=int((d.LA_high_RMSE<d.DS_high_RMSE).sum()),
        additional_spectra=dict(n=len(ext),LA_lower_regional_RMSE=int((ext.LA_contiguous_high_RMSE<ext.Voigt_contiguous_high_RMSE).sum()),
            median_RMSE_reduction_percent=float(ext.contiguous_improvement_percent.median()),
            minimum_reduction_percent=float(ext.contiguous_improvement_percent.min()),
            maximum_reduction_percent=float(ext.contiguous_improvement_percent.max())),
        synthetic=dict(n=len(paired),free_lower_RMSE=int(lower.sum()),free_larger_absolute_error=int(larger.sum()),
            lower_RMSE_and_larger_error=int((lower&larger).sum()),
            max_MCMC_absolute_error_pp=float(abs(me).max()),error_by_panel_and_noise=records),
        joint_background=dict(pairs=len(j),accepted=len(accepted),unsuccessful=int((~ok).sum()),
            accepted_with_anchor_boundary=int(boundary.sum()),
            accepted_full_RMSE_lower=int((accepted.RMSE_joint<accepted.RMSE_fixed).sum()),
            accepted_high_RMSE_lower=int((accepted.high_RMSE_joint<accepted.high_RMSE_fixed).sum()),
            accepted_high_RMSE_higher=int((accepted.high_RMSE_joint>accepted.high_RMSE_fixed).sum())))

if __name__=='__main__':
    out=ROOT/'results';out.mkdir(exist_ok=True)
    result=calculate()
    (out/'publication_statistics.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
