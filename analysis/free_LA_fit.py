"""Independent ten-peak LA/Voigt least squares on a persisted fixed background.

Only Ce(IV) v/u use LA(alpha=.9,beta=2); each peak has independent center,
area, sigma and gamma. Profiles have unit area in the observed window.  There
are no spin-orbit, multiplet-area, or shared-width parameter relationships.
The extra constrained-LS seed is data-derived, never a ground-truth seed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
import numpy as np
import pandas as pd
from scipy.fft import irfft,rfft,rfftfreq
from scipy.optimize import least_squares
from scipy.special import wofz

SOURCE=Path(__file__).resolve().parents[1]
WORK=SOURCE/'results/synthetic_refits'
DATA=WORK
sys.path.insert(0,str(SOURCE))
from ce3d_xps_model import Ce3dModel,MAIN_CENTERS,MAIN_WEIGHTS,GAMMA_INDICES,MAIN_COMPONENT_NAMES,LA_GRID_N,LA_GRID_DX

def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');os.replace(tmp,path)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class FreeLAModel:
    def __init__(self,e,area_scale):
        self.e=np.asarray(e,float);self.area_scale=float(area_scale)
        self.grid=(np.arange(LA_GRID_N)-LA_GRID_N//2)*LA_GRID_DX
        self.freq=rfftfreq(LA_GRID_N,LA_GRID_DX)
        self.power=np.where(self.grid>=0,.9,2.)
        self.last=None
    def evaluate(self,q):
        q=np.asarray(q,float)
        if self.last is not None and np.array_equal(q,self.last[0]):return self.last[1:]
        t=q.reshape(10,4);curves=[];jac=np.zeros((len(self.e),40))
        for j,(offset,area,sigma,gamma) in enumerate(t):
            z=self.e-MAIN_CENTERS[j]-offset
            if j in (4,5):
                rr=(self.grid/gamma)**2
                h=(1+rr)**(-self.power)
                dh=h*(2*self.power*rr/(gamma*(1+rr)))
                H=rfft(h);G=np.exp(-2*np.pi**2*sigma**2*self.freq**2)
                conv=irfft(H*G,n=LA_GRID_N)
                conv_sigma=irfft(H*G*(-4*np.pi**2*sigma*self.freq**2),n=LA_GRID_N)
                conv_gamma=irfft(rfft(dh)*G,n=LA_GRID_N)
                p=np.interp(z,self.grid,conv)
                idx=np.clip(np.searchsorted(self.grid,z,side='right')-1,0,LA_GRID_N-2)
                dc=-(conv[idx+1]-conv[idx])/LA_GRID_DX
                dc[(z<self.grid[0])|(z>self.grid[-1])]=0
                ds=np.interp(z,self.grid,conv_sigma);dg=np.interp(z,self.grid,conv_gamma)
            else:
                zz=(z+1j*gamma)/(sigma*np.sqrt(2))
                w=wofz(zz);wp=-2*zz*w+2j/np.sqrt(np.pi)
                den=sigma*np.sqrt(2*np.pi)
                p=w.real/den
                dc=(-wp/(sigma*np.sqrt(2))).real/den
                ds=-p/sigma+(-zz*wp/sigma).real/den
                dg=(1j*wp/(sigma*np.sqrt(2))).real/den
            norm=np.trapz(p,self.e)
            pn=p/norm
            deriv=np.column_stack((dc,ds,dg))
            deriv=(deriv-pn[:,None]*np.trapz(deriv,self.e,axis=0)[None,:])/norm
            amp=area*self.area_scale
            curves.append(amp*pn)
            jac[:,4*j]=amp*deriv[:,0];jac[:,4*j+1]=self.area_scale*pn
            jac[:,4*j+2]=amp*deriv[:,1];jac[:,4*j+3]=amp*deriv[:,2]
        comps=np.column_stack(curves);prediction=comps.sum(axis=1)
        self.last=(q.copy(),prediction,jac,comps)
        return prediction,jac,comps
    def predict(self,q):return self.evaluate(q)[0]
    def jacobian(self,q):return self.evaluate(q)[1]

def expand_constrained(t,area_scale):
    amps=t[3]*MAIN_WEIGHTS*np.where(np.arange(10)<4,t[2],1-t[2])
    return np.column_stack((np.full(10,t[0]),amps/area_scale,np.full(10,t[1]),t[GAMMA_INDICES])).ravel()

def check_model():
    e=np.round(np.arange(870.58,925,.1),2);ref=Ce3dModel(e,np.zeros_like(e),'LA')
    t=np.array([.123,.81,.31,13000.,.43,.61,.22,.92,.31,0.,0.])
    free=FreeLAModel(e,13000.);q=expand_constrained(t,13000.)
    actual=free.predict(q);expected=ref.predict(t)[0]
    np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=1e-11)
    jac=free.jacobian(q).copy();numeric=np.empty_like(jac)
    for j in range(40):
        h=1e-6;qp=q.copy();qm=q.copy();qp[j]+=h;qm[j]-=h
        numeric[:,j]=(free.predict(qp)-free.predict(qm))/(2*h)
    errs=np.max(abs(jac-numeric),axis=0)/np.maximum(np.max(abs(numeric),axis=0),1.)
    assert errs.max()<2e-6,errs
    out=dict(status='passed',prediction_max_abs_difference=float(np.max(abs(actual-expected))),
        jacobian_max_relative_column_error=float(errs.max()),points=len(e),
        normalization='observed window',free_parameters=40,
        lineshape='LA(alpha=.9,beta=2) only at component indices 4/5, eight other Voigt')
    dump(WORK/'free_LA_model_check.json',out);return out

def fit_case(case_id):
    manifest=json.loads((WORK/'benchmark_manifest.json').read_text(encoding='utf-8'))
    c=next(c for c in manifest['cases'] if c['case_id']==case_id)
    out=DATA/'free_LA'/case_id;out.mkdir(parents=True,exist_ok=True)
    if (out/'fit_result.json').exists():
        r=json.loads((out/'fit_result.json').read_text(encoding='utf-8'))
        if r.get('status')=='completed':return r
    start=time.time();d=pd.read_csv(c['input']);e=d.energy_eV.to_numpy(float)
    raw=d.raw.to_numpy(float);bg=d.background.to_numpy(float);net=raw-bg
    area_scale=max(float(np.trapz(np.maximum(net,0),e)),100.)
    model=FreeLAModel(e,area_scale)
    ref=Ce3dModel(e,net,'LA');t=ref.initial();lo,hi=ref.bounds
    # Common ordinary-SSE constrained control; no prior or noise weighting.
    rscale=max(float(np.max(abs(net))),1.)
    ctl=least_squares(lambda v:(ref.predict(np.r_[v,t[9:]])[0]-net)/rscale,
        np.clip(t[:9],lo[:9]+1e-8,hi[:9]-1e-8),bounds=(lo[:9]+1e-8,hi[:9]-1e-8),
        x_scale=ref.scales[:9],max_nfev=450,ftol=1e-10,xtol=1e-10,gtol=1e-10)
    tc=np.r_[ctl.x,t[9:]];ctl_pred=ref.predict(tc)[0]
    control=dict(parameters=dict(zip(ref.param_names,tc.tolist())),SSE=float(np.sum((net-ctl_pred)**2)),
        RMSE=float(np.sqrt(np.mean((net-ctl_pred)**2))),Ce3_percent=float(100*tc[2]),success=bool(ctl.success),nfev=int(ctl.nfev),
        method='constrained LA ordinary least squares; diagnostic seed only; no ground-truth initialization')
    dump(out/'constrained_LS_control.json',control)
    q0=expand_constrained(t,area_scale)
    boundslo=np.tile([-1.,0.,.05,.1],10)
    boundshi=np.tile([1.,5.,1.5,2.],10)
    seed=int.from_bytes(hashlib.sha256(case_id.encode('utf-8')).digest()[:4],'little')
    rng=np.random.default_rng(seed)
    starts=[('default_data_derived',q0),('constrained_LS_seed',expand_constrained(tc,area_scale))]
    for k in range(2):
        qp=q0.reshape(10,4).copy()
        qp[:,0]=rng.uniform(-.5,.5,10);qp[:,1]*=rng.uniform(.5,1.7,10)
        qp[:,2]=rng.uniform(.35,1.4,10);qp[:,3]=rng.uniform(.15,1.1,10)
        starts.append((f'fixed_random_start_{k+1}',qp.ravel()))
    records=[];best=None
    for k,(label,initial) in enumerate(starts):
        initial=np.clip(initial,boundslo+1e-10,boundshi-1e-10)
        opt=least_squares(lambda v:(model.predict(v)-net)/rscale,initial,
            jac=lambda v:model.jacobian(v)/rscale,bounds=(boundslo,boundshi),
            x_scale=np.tile([.3,.1,.5,.5],10),loss='linear',max_nfev=1000,
            ftol=1e-9,xtol=1e-9,gtol=1e-9)
        pred=model.predict(opt.x);sse=float(np.sum((net-pred)**2))
        records.append(dict(start=k,label=label,SSE=sse,success=bool(opt.success),status=int(opt.status),
            nfev=int(opt.nfev),njev=int(opt.njev),optimality=float(opt.optimality),message=str(opt.message),parameters=opt.x.tolist()))
        if best is None or sse<best[0]:best=(sse,opt.x.copy(),k)
        dump(out/'start_results.json',records)
        print(case_id,'start',k,'SSE',sse,'nfev',opt.nfev,flush=True)
    sse,q,k=best;pred,jac,comps=model.evaluate(q);pp=q.reshape(10,4)
    areas=np.trapz(comps,e,axis=0);x=areas[:4].sum()/areas.sum()
    assert sse<=control['SSE']*(1+1e-7)+1e-6,(sse,control['SSE'])
    ptable=pd.DataFrame(dict(component=MAIN_COMPONENT_NAMES,center_eV=MAIN_CENTERS+pp[:,0],
        window_area=areas,sigma_eV=pp[:,2],gamma_parameter_eV=pp[:,3],
        gamma_definition=['LA kernel scale' if j in (4,5) else 'Lorentzian HWHM' for j in range(10)]))
    ptable.to_csv(out/'peak_table.csv',index=False)
    curves=pd.DataFrame(dict(energy_eV=e,raw=raw,background=bg,net=net,net_fit=pred,total_fit=bg+pred,residual=net-pred))
    for j,name in enumerate(MAIN_COMPONENT_NAMES):curves[name]=comps[:,j]
    curves.to_csv(out/'fit_curves.csv',index=False)
    result=dict(case_id=case_id,panel=c['panel'],status='completed',method='independent LA/Voigt ordinary least squares',
        free_parameter_count=40,Ce3_percent=float(100*x),Ce3_percent_true=c['Ce3_percent_true_nominal'],
        sigma_true_eV=c['sigma_true_eV'],noise_percent=c['noise_percent'],
        RMSE=float(np.sqrt(sse/len(e))),SSE=sse,best_start=k,best_start_success=records[k]['success'],
        starts=4,random_seed=seed,constrained_LS_control=control,
        residual_objective='unweighted ordinary SSE; loss=linear',
        area_definition='CeIII component areas / all component areas within identical observed window',
        position_bounds='each nominal labeled center +/-1.0 eV',width_bounds=dict(sigma_eV=[.05,1.5],gamma_eV=[.1,2.]),
        parameter_relationships='none; independent center, area, sigma and gamma for all ten labeled peaks',
        LA_alpha=.9,LA_beta=2.,asymmetric_components=[MAIN_COMPONENT_NAMES[4],MAIN_COMPONENT_NAMES[5]],
        prepared_input=str(c['input']),prepared_input_sha256=sha(c['input']),
        background_sha256=hashlib.sha256(np.ascontiguousarray(bg).tobytes()).hexdigest(),
        elapsed_seconds=time.time()-start,run_dir=str(out),posterior_interval=None)
    dump(out/'fit_result.json',result)
    return result

def collect():
    rows=[]
    for p in (DATA/'free_LA').glob('*/fit_result.json'):
        r=json.loads(p.read_text(encoding='utf-8'))
        if r.get('status')!='completed':continue
        rows.append({k:r[k] for k in ['case_id','panel','Ce3_percent','Ce3_percent_true','sigma_true_eV','noise_percent','RMSE','SSE','best_start','best_start_success','elapsed_seconds','run_dir','prepared_input_sha256','background_sha256']}
            | dict(constrained_LS_RMSE=r['constrained_LS_control']['RMSE'],constrained_LS_Ce3_percent=r['constrained_LS_control']['Ce3_percent']))
    if rows:
        df=pd.DataFrame(rows).sort_values('case_id');tmp=WORK/f'free_LA_summary_{os.getpid()}.tmp.csv'
        df.to_csv(tmp,index=False);os.replace(tmp,WORK/'free_LA_summary.csv')
    return rows

def batch(workers):
    manifest=json.loads((WORK/'benchmark_manifest.json').read_text(encoding='utf-8'))
    queue=[c['case_id'] for c in manifest['cases'] if not (DATA/'free_LA'/c['case_id']/'fit_result.json').exists()]
    active={};failed=[];start=time.time();env=os.environ.copy()
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONIOENCODING='utf-8')
    logs=DATA/'free_LA_logs';logs.mkdir(parents=True,exist_ok=True)
    while queue or active:
        while queue and len(active)<workers:
            key=queue.pop(0);f=(logs/f'{key}.log').open('a',encoding='utf-8')
            p=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--case',key],stdout=f,stderr=subprocess.STDOUT,
                env=env,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));active[key]=(p,f)
        for key,(p,f) in list(active.items()):
            code=p.poll()
            if code is not None:
                f.close();del active[key]
                if code:failed.append(dict(case_id=key,exit=code))
                print('FINISHED',key,'exit',code,flush=True);collect()
        dump(WORK/'free_LA_batch_status.json',dict(running={k:p.pid for k,(p,f) in active.items()},pending=len(queue),completed=len(collect()),failures=failed,elapsed_seconds=time.time()-start))
        if queue or active:time.sleep(5)
    print('FREE LA COMPLETE',len(collect()),'failed',failed,flush=True)

if __name__=='__main__':
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');ap.add_argument('--case');ap.add_argument('--batch',action='store_true');ap.add_argument('--workers',type=int,default=2);ap.add_argument('--collect',action='store_true');args=ap.parse_args()
    if args.check:print(json.dumps(check_model(),indent=2))
    elif args.case:
        try:print(json.dumps(fit_case(args.case),ensure_ascii=False,indent=2))
        except Exception:traceback.print_exc();sys.exit(1)
    elif args.batch:batch(args.workers)
    elif args.collect:print('COLLECT',len(collect()))
