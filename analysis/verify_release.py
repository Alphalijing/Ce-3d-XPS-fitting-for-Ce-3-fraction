"""Read-only checks of the distributed evidence; writes a report under results/."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pandas as pd
from summarize_results import calculate

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ce3d_io import load_spectrum

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify():
    manifest=ROOT/'MANIFEST.json'
    hashes=0
    if manifest.exists():
        for r in json.loads(manifest.read_text(encoding='utf-8'))['files']:
            p=ROOT/r['file'];assert p.is_file(),r['file']
            assert p.stat().st_size==r['bytes'],r['file']
            assert sha(p)==r['sha256'],r['file']
            hashes+=1
    idx=pd.read_csv(ROOT/'data/experimental/la_series/sample_index.csv')
    assert idx.La_at_percent.tolist()==[0,1,2,4,10,20]
    main=pd.read_csv(ROOT/'data/experimental/la_series/Figure4_all_model_results.csv')
    assert len(main)==18
    for r in idx.itertuples():
        src=ROOT/r.raw_file;assert sha(src)==r.raw_file_sha256
        spec=load_spectrum(src)
        curves=pd.read_csv(ROOT/f'data/experimental/la_series/curves/{r.sample_id}.csv')
        np.testing.assert_allclose(spec.raw,curves.raw,rtol=0,atol=1e-7)
        np.testing.assert_allclose(spec.background,curves.LA_fixed_background,rtol=0,atol=1e-7)
        for model in ('LA','Voigt','DS'):
            values=curves.loc[curves.high_mask,'raw']-curves.loc[curves.high_mask,model+'_fixed_total']
            rmse=float(np.sqrt(np.mean(values**2)))
            reference=main[(main.La_at_percent==r.La_at_percent)&(main['shape']==model)].iloc[0]
            np.testing.assert_allclose(rmse,reference.high_RMSE,rtol=1e-10)
    posterior=json.loads((ROOT/'data/experimental/la_series/posterior/index.json').read_text())
    assert len(posterior)==12
    for r in posterior:
        with np.load(ROOT/r['file']) as f:
            draws=f['Ce3_fraction'];assert len(draws)==320000
            q=np.quantile(draws*100,[.025,.5,.975])
        row=main[(np.isclose(main['sample'],float(r['sample'])))&(main['shape']==r['shape'])].iloc[0]
        np.testing.assert_allclose(q,[row.Ce3_low,row.Ce3_window_percent,row.Ce3_high],rtol=0,atol=1e-8)
    ext=pd.read_csv(ROOT/'data/experimental/additional_19/summary.csv')
    expected={f'batch2_{n}' for n in [1,3,4,5,6,7,*range(9,22)]}
    assert set(ext['sample'])==expected
    for r in ext.itertuples():
        assert sha(ROOT/r.source)==r.raw_file_sha256
        c=pd.read_csv(ROOT/f'data/experimental/additional_19/curves/{r.sample}.csv')
        fixed=pd.read_csv(ROOT/r.prepared_file)
        np.testing.assert_allclose(fixed.background,c.noncrossing_Shirley,rtol=0,atol=1e-8)
        assert int(c.high_mask.sum())==int(r.contiguous_high_points)
        for shape in ('LA','Voigt'):
            rmse=np.sqrt(np.mean((c.loc[c.high_mask,'raw']-c.loc[c.high_mask,shape+'_total'])**2))
            np.testing.assert_allclose(rmse,getattr(r,shape+'_contiguous_high_RMSE'),rtol=1e-9)
    bench=json.loads((ROOT/'data/synthetic/benchmark_manifest.json').read_text())
    assert len(bench['cases'])==48
    for case in bench['cases']:
        p=ROOT/case['input'];assert sha(p)==case['prepared_sha256']
        raw=ROOT/case['original'];assert sha(raw)==case['original_sha256']
        inp=pd.read_csv(p);truth=pd.read_csv(ROOT/case['clean_reference'])
        t=truth[truth.in_fit_window==1]
        assert len(inp)==545
        np.testing.assert_allclose(inp.background,t.true_background,rtol=0,atol=1e-9)
        np.testing.assert_array_equal(inp.raw,t.raw_noisy)
        if case['panel']=='a':assert np.all(inp.background==0)
    st=calculate()
    assert st['additional_spectra']['LA_lower_regional_RMSE']==19
    assert st['synthetic']['lower_RMSE_and_larger_error']==47
    assert st['joint_background']['accepted']==17
    assert st['joint_background']['accepted_with_anchor_boundary']==15
    # Explicit import of SVG backend catches the earlier frozen-build omission.
    import matplotlib.backends.backend_svg
    return dict(status='passed',files_hash_checked=hashes,main_samples=6,
        full_Ce3_posteriors_verified=12,retained_Ce3_draws_verified=3840000,
        additional_spectra=19,synthetic_cases=48,MCMC_started=False,statistics=st)

if __name__=='__main__':
    result=verify();out=ROOT/'results';out.mkdir(exist_ok=True)
    (out/'release_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='statistics'},indent=2))
