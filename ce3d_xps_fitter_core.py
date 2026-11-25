"""Ce 3d XPS Fitter 2.0: fixed background, constrained peaks, frozen-proposal MCMC."""
from __future__ import annotations
import argparse
from datetime import datetime
import json, os, shutil, time, traceback, sys
from pathlib import Path
from threading import Event
from typing import Callable
import emcee
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize
from scipy.special import logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends import backend_svg as _backend_svg  # Required export format; validate before sampling.
from ce3d_io import load_spectrum
from ce3d_vamas import json_value, source_record
from ce3d_xps_model import Ce3dModel
from ce3d_results import write_result_report, organize_outputs, write_batch_report

__version__='2.0.1'
DEFAULT_WALKERS=80
DEFAULT_STEPS=8000
DEFAULT_BURN=4000
CHECK_EVERY=2000
PILOT_STEPS=2000

def write_json(path, value):
    path=Path(path); temp=path.with_suffix('.tmp.json')
    temp.write_text(json.dumps(json_value(value),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    os.replace(temp,path)

def _emit(log, message):
    if log: log(message)

def _progress(callback, **event):
    if callback: callback(event)

class CancelledBeforeSampling(Exception):pass

class ChainStore:
    """Only closed, atomically published chunks are ever read."""
    def __init__(self, folder):
        self.folder=Path(folder)/'chain';self.folder.mkdir();self.blocks=[];self.steps=0
    def append(self, chain, logp, state):
        if len(chain)==0:return
        end=self.steps+len(chain)
        path=self.folder/f'chain_{self.steps:07d}_{end:07d}.npz'
        temp=path.with_suffix('.tmp')
        with temp.open('wb') as f:np.savez(f,chain=chain,log_prob=logp)
        os.replace(temp,path)
        self.blocks.append(dict(start=self.steps,end=end,file=path.name));self.steps=end
        write_json(self.folder/'manifest.json',dict(steps=end,blocks=self.blocks))
        state_path=self.folder/'last_state.npz';temp=state_path.with_suffix('.tmp')
        with temp.open('wb') as f:np.savez(f,coords=state.coords,log_prob=state.log_prob)
        os.replace(temp,state_path)
        rs=state.random_state
        if rs is not None:write_json(self.folder/'random_state.json',dict(algorithm=rs[0],keys=rs[1],position=rs[2],has_gauss=rs[3],cached_gaussian=rs[4]))
    def read(self):
        chains=[];lps=[]
        for b in self.blocks:
            with np.load(self.folder/b['file'],allow_pickle=False) as z:
                chains.append(z['chain']);lps.append(z['log_prob'])
        return np.concatenate(chains),np.concatenate(lps)

def initialize(model, seed, stop_event=None):
    """Fresh-data least-squares seed, then four starts of the unchanged posterior."""
    rng=np.random.default_rng(seed);initial=model.initial().copy()
    lo,hi=(np.array(v,dtype=float) for v in model.bounds)
    eps=np.maximum((hi-lo)*1e-9,1e-9)
    indices=np.r_[np.arange(9),np.arange(11,model.ndim)].astype(int)
    signal_scale=max(float(np.max(abs(model.net))),1.)
    def seed_residual(v):
        t=initial.copy();t[indices]=v
        if not np.isfinite(model.log_prior(t)):return np.full(len(model.energy),1e6)
        return (model.predict(t)[0]-model.net)/signal_scale
    ls_candidates=[]
    for j in range(3):
        if stop_event is not None and stop_event.is_set():raise CancelledBeforeSampling()
        t=initial.copy()
        if j:
            t[0]=rng.uniform(-.3,.3);t[1]=rng.uniform(.35,.9);t[2]=rng.uniform(.08,.35);t[4:9]=rng.uniform(.2,1.1,5)
        opt=least_squares(seed_residual,np.clip(t[indices],lo[indices]+eps[indices],hi[indices]-eps[indices]),
            bounds=(lo[indices]+eps[indices],hi[indices]-eps[indices]),x_scale=model.scales[indices],max_nfev=450)
        candidate=initial.copy();candidate[indices]=opt.x
        if np.isfinite(model.log_prior(candidate)):ls_candidates.append((float(opt.fun@opt.fun),candidate))
    if ls_candidates:initial=min(ls_candidates,key=lambda v:v[0])[1]
    widths=np.r_[.15,.12,.05,max(initial[3]*.06,1.),[.12]*5,.4,.4]
    if model.satellites:widths=np.r_[widths,.006,.006,.04,.025]
    candidates=[]
    for k,(a,b) in enumerate([(1e4,300.),(1e6,.01),(5e5,30.),(1e4,1.)]):
        if stop_event is not None and stop_event.is_set():raise CancelledBeforeSampling()
        q=initial.copy()+rng.normal(size=model.ndim)*widths
        q[9:11]=np.log([a,b])+rng.normal(size=2)*.2
        q=np.clip(q,lo+eps,hi-eps)
        if not np.isfinite(model.log_prior(q)):q=initial.copy();q[9:11]=np.log([a,b])
        def objective(v):
            value=model.log_prob(v*model.scales)[0]
            return -float(value) if np.isfinite(value) else 1e100
        opt=minimize(objective,q/model.scales,method='L-BFGS-B',
            bounds=list(zip((lo+eps)/model.scales,(hi-eps)/model.scales)),
            options=dict(maxiter=350,ftol=1e-10,maxls=40))
        qopt=opt.x*model.scales
        if not np.isfinite(model.log_prob(qopt)[0]):qopt=q
        candidates.append(dict(start=k,parameters=qopt.tolist(),log_posterior=float(model.log_prob(qopt)[0]),
            success=bool(opt.success),message=str(opt.message)))
    best=np.array(max(candidates,key=lambda c:c['log_posterior'])['parameters'])
    alternatives=[c for c in candidates if c['parameters'][10]<0]
    alt=np.array(max(alternatives,key=lambda c:c['log_posterior'])['parameters']) if alternatives else initial.copy()
    if not alternatives:alt[9:11]=np.log([1e6,.01])
    widths=np.r_[.0015,.025,.015,max(initial[3]*.03,.1),.06,.06,.035,.06,.035,.5,.5]
    if model.satellites:widths=np.r_[widths,.004,.004,.03,.015]
    return rng,best,alt,candidates,widths,dict(least_squares_starts=len(ls_candidates),least_squares_seed=initial.tolist())

def initialize_walkers(model,rng,best,alt,candidates,widths,walkers):
    result=[]
    for i in range(walkers):
        center=best if i<walkers//2 else alt if i<3*walkers//4 else np.asarray(candidates[i%len(candidates)]['parameters'])
        for _ in range(100000):
            q=center+rng.normal(size=model.ndim)*widths
            if np.isfinite(model.log_prob(q)[0]):break
        else:raise RuntimeError('Could not construct finite-prior initial walkers.')
        result.append(q)
    return np.array(result)

def freeze_proposal(pilot, model, best, alt, widths):
    train=pilot[len(pilot)//2:].reshape(-1,model.ndim);props=[];counts=[]
    for k,mask in enumerate([train[:,10]>=0,train[:,10]<0]):
        y=train[mask]/model.scales;counts.append(len(y))
        if len(y)>=2*model.ndim:
            mean=y.mean(axis=0);cov=np.cov(y,rowvar=False)
            cov=2*(.95*cov+.05*np.diag(np.diag(cov))+np.eye(model.ndim)*1e-10)
        else:
            mean=(best if k==0 else alt)/model.scales
            cov=np.diag(np.maximum(widths/model.scales,.01)**2)*8
        chol=np.linalg.cholesky(cov)
        props.append((mean,chol,np.linalg.inv(chol),float(np.log(np.diag(chol)).sum())))
    weights=np.array([.75,.25])
    def logq(coords):
        y=coords/model.scales;terms=[]
        for w,(mean,chol,inv,ld) in zip(weights,props):
            z=(y-mean)@inv.T
            terms.append(np.log(w)-.5*model.ndim*np.log(2*np.pi)-ld-.5*np.sum(z*z,axis=1))
        return logsumexp(terms,axis=0)
    def proposal(coords,random):
        choices=random.choice(2,len(coords),p=weights);z=random.normal(size=coords.shape);new=np.empty_like(coords)
        for k,(mean,chol,_,_) in enumerate(props):
            mask=choices==k;new[mask]=(mean+z[mask]@chol.T)*model.scales
        return new,logq(coords)-logq(new)
    return proposal,dict(frozen=True,weights=weights,training_counts=counts,
        means_scaled=[v[0] for v in props],covariances_scaled=[v[1]@v[1].T for v in props],
        full_mixture_Hastings_correction=True,untruncated_normal=True,scale=model.scales)

def convergence_diagnostics(chain,lps,burn,history):
    n=len(chain);post=chain[burn:]
    if len(post)<2:return dict(evaluated=False,passed=False,reason='burn-in not complete')
    tau=np.asarray(emcee.autocorr.integrated_time(post,tol=0,quiet=True))
    tlp=float(emcee.autocorr.integrated_time(lps[burn:,:,None],tol=0,quiet=True)[0])
    all_tau=np.r_[tau,tlp]
    recent=history[-2:]+[dict(all_tau=all_tau.tolist())]
    changes=[float(np.max(abs(np.asarray(b['all_tau'])/np.asarray(a['all_tau'])-1))) for a,b in zip(recent[:-1],recent[1:])]
    windows=[]
    for a in range(max(burn,n-6000),n,2000):
        if n-a<2000:continue
        windows.append(dict(start=a,end=a+2000,Ce3_quantiles_percent=(100*np.quantile(chain[a:a+2000,:,2],[.025,.5,.975])).tolist()))
    ranges=np.ptp([w['Ce3_quantiles_percent'] for w in windows],axis=0) if len(windows)==3 else None
    finite=bool(np.isfinite(all_tau).all() and (all_tau>0).all())
    gates=dict(positive_finite_tau=finite,
        retained_length_at_least_100_max_tau=bool(finite and len(post)>=100*max(all_tau)),
        three_checks_all_tau_change_below_5_percent=bool(len(recent)==3 and np.isfinite(changes).all() and max(changes)<.05),
        three_window_quantile_range_below_0p1_pp=bool(ranges is not None and max(ranges)<.1),
        minimum_total_steps_10000=bool(n>=10000))
    return dict(evaluated=True,passed=bool(all(gates.values())),steps=n,retained_steps=len(post),
        tau_parameters=tau.tolist(),tau_log_posterior=tlp,all_tau=all_tau.tolist(),
        max_tau=float(max(all_tau)) if finite else None,
        retained_over_max_tau=float(len(post)/max(all_tau)) if finite else None,
        tau_changes_max_relative=changes,Ce3_windows=windows,window_quantile_ranges_pp=ranges,
        checks=gates,low_b_sample_fraction=float(np.mean(post[:,:,10]<0)),
        low_b_walkers_last_2000=int((np.median(chain[-min(n,2000):,:,10],axis=0)<0).sum()),
        dE_fraction_below_minus_0p99=float(np.mean(post[:,:,0]<-.99)))

def _save_background(out,spectrum,label='Fixed Shirley background'):
    e,y,bg=spectrum.energy,spectrum.raw,spectrum.background
    pd.DataFrame(dict(energy_eV=e,raw=y,background=bg,net=y-bg)).to_csv(out/'prepared_spectrum.csv',index=False)
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    ax.plot(e,y,color='black',lw=.9,label='Measured spectrum');ax.plot(e,bg,label=label,color='#009E73')
    ax.invert_xaxis();ax.set_xlabel('Binding energy (eV)');ax.set_ylabel('Intensity (input units)');ax.legend()
    fig.savefig(out/'background.png',dpi=180);plt.close(fig)

def _export(out,spectrum,model,store,burn,status,config,diagnostics,elapsed):
    result=dict(software='Ce 3d XPS Fitter',version=__version__,status=status,configuration=config,
        total_steps=store.steps,burn=burn,elapsed_seconds=elapsed,
        experiment_metadata_file='experiment_metadata.json',convergence=diagnostics,
        model=model.metadata(),uncertainty_scope='Conditional on the selected fixed background, peaks, priors and noise likelihood; does not include all model or experimental uncertainty.')
    if store.steps<=burn:
        result['posterior_available']=False;write_json(out/'fit_result.json',result)
        write_result_report(out,result)
        return result
    chain,lps=store.read();post=chain[burn:];flat=post.reshape(-1,model.ndim)
    q=100*np.quantile(flat[:,2],[.025,.5,.975]);theta=np.median(flat,axis=0)
    comps=model.components(theta);netfit=comps.sum(axis=1);total=spectrum.background+netfit;resid=spectrum.raw-total
    pd.DataFrame(flat,columns=model.param_names).to_csv(out/'posterior_samples_ce3d.csv',index=False)
    quant=np.quantile(flat,[.025,.5,.975],axis=0).T
    pd.DataFrame(dict(parameter=model.param_names,q025=quant[:,0],median=quant[:,1],q975=quant[:,2],
        mean=flat.mean(axis=0),std=flat.std(axis=0))).to_csv(out/'posterior_parameters.csv',index=False)
    pd.DataFrame(model.peak_table(theta)).to_csv(out/'peak_table.csv',index=False)
    curves=dict(energy_eV=spectrum.energy,raw=spectrum.raw,background=spectrum.background,
        net=spectrum.raw-spectrum.background,total_fit=total,net_fit=netfit,residual=resid)
    for j,name in enumerate(model.component_names):curves[name]=comps[:,j]
    pd.DataFrame(curves).to_csv(out/'fit_curves.csv',index=False)
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True,layout='constrained',gridspec_kw={'height_ratios':[3,1]})
    e=spectrum.energy;bg=spectrum.background
    axes[0].plot(e,spectrum.raw,'k',lw=.8,label='Measured spectrum');axes[0].plot(e,total,color='#0072B2',label='Total fit')
    axes[0].plot(e,bg,color='#009E73',ls='--',label='Fixed background')
    for j,name in enumerate(model.component_names):
        color='#009E73' if j<4 else '#CC79A7' if j>=10 else '#D55E00'
        label='Ce(III)' if j==0 else 'Ce(IV)' if j==4 else 'Ce(IV) satellites' if j==10 else None
        axes[0].plot(e,bg+comps[:,j],color=color,lw=.6,alpha=.8,label=label)
        axes[0].fill_between(e,bg,bg+comps[:,j],color=color,alpha=.07)
    background_label='fixed Shirley background' if config['background']=='fixed_shirley' else 'provided background'
    axes[0].set_title(f'{model.lineshape} + {background_label} | Ce(III) = {q[1]:.2f}% [{q[0]:.2f}, {q[2]:.2f}]%')
    axes[0].set_ylabel('Intensity (input units)');axes[0].legend(fontsize=8)
    axes[1].plot(e,resid,'k',lw=.7);axes[1].axhline(0,color='gray',lw=.6)
    axes[1].set_xlabel('Binding energy (eV)');axes[1].set_ylabel('Residual');axes[1].invert_xaxis()
    fig.savefig(out/'fit_ce3d.png',dpi=180);fig.savefig(out/'fit_ce3d.svg');plt.close(fig)
    acceptance=np.any(post[1:]!=post[:-1],axis=2).mean(axis=0) if len(post)>1 else np.full(chain.shape[1],np.nan)
    result.update(posterior_available=True,Ce3_fraction_median=float(np.median(flat[:,2])),
        Ce3_fraction_mean=float(np.mean(flat[:,2])),Ce3_fraction_std=float(np.std(flat[:,2])),
        Ce3_percent=dict(q025=float(q[0]),median=float(q[1]),q975=float(q[2])),
        posterior_median=dict(zip(model.param_names,theta)),
        Ce3_area_definition='Ce(III) component areas / all Ce component areas within fitted window; Ce(IV) satellites included when enabled',
        metrics=dict(RMSE=float(np.sqrt(np.mean(resid**2))),SSE=float(resid@resid)),
        production_acceptance_median=float(np.nanmedian(acceptance)),
        bounds_diagnostics=dict(dE_fraction_below_minus_0p99=float(np.mean(flat[:,0]<-.99))))
    write_json(out/'fit_result.json',result)
    write_result_report(out,result,flat=flat,param_names=model.param_names)
    return result

def _fit_one(path,outroot,*,model, satellites,mode,walkers,steps,burn,seed,background_mode,
             vamas_block,vamas_ordinate,kinetic_offset_ev,satellite_gamma,max_steps,log,progress,stop_event,
             file_index,file_count):
    out=outroot/f"{Path(path).stem}_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}";out.mkdir(parents=True)
    start=time.time();store=None
    _emit(log,f'Starting: {Path(path).name} | {model} + fixed Shirley background | '+('satellites included' if satellites else 'satellites excluded')+' | '+('Quick fit' if mode=='quick' else 'Convergence validation'))
    _progress(progress,event='start',steps=0,total_steps=steps if mode=='quick' else None,elapsed_seconds=0,
              file_index=file_index,file_count=file_count,output_dir=str(out))
    config=dict(version=__version__,source=str(Path(path).resolve()),model=model,satellites=bool(satellites),mode=mode,
        walkers=walkers,steps=steps if mode=='quick' else None,burn=burn,pilot_steps=PILOT_STEPS,
        checkpoint_every=CHECK_EVERY,seed=seed,background=background_mode,satellite_gamma_HWHM_eV=satellite_gamma,
        max_steps=max_steps,proposal_moves=dict(DEMove=.5,DESnookerMove=.2,frozen_independence_MH=.3),
        LA_alpha=.9 if model=='LA' else None,LA_beta=2. if model=='LA' else None)
    write_json(out/'run_configuration.json',config)
    write_json(out/'experiment_metadata.json',dict(schema_version='2.0',source=source_record(path),
        experiment_conditions={},status='awaiting_input_parsing'))
    try:
        spectrum=load_spectrum(path,background_mode=background_mode,vamas_block=vamas_block,
            vamas_ordinate=vamas_ordinate,kinetic_offset_ev=kinetic_offset_ev)
        write_json(out/'experiment_metadata.json',spectrum.metadata)
        _save_background(out,spectrum,'Fixed Shirley background' if background_mode=='fixed_shirley' else 'Provided background')
        for note in spectrum.notes:_emit(log,str(note))
        target=Ce3dModel(spectrum.energy,spectrum.raw-spectrum.background,lineshape=model,satellites=satellites,satellite_gamma=satellite_gamma)
        if walkers<2*target.ndim:raise ValueError(f'walkers must be >= {2*target.ndim} for this model.')
        _emit(log,'Initializing peak parameters and sampler; background plot and acquisition metadata JSON saved.')
        init_start=time.time();rng,best,alt,candidates,widths,seed_info=initialize(target,seed,stop_event)
        write_json(out/'initialization.json',dict(optimizations=candidates,seconds=time.time()-init_start,**seed_info))
        if stop_event is not None and stop_event.is_set():raise CancelledBeforeSampling()
        walk=initialize_walkers(target,rng,best,alt,candidates,widths,walkers)
        np.save(out/'initial_positions.npy',walk)
        np.random.seed(seed);store=ChainStore(out)
        sampler=emcee.EnsembleSampler(walkers,target.ndim,target.log_prob,vectorize=True,
            moves=[(emcee.moves.DEMove(),.8),(emcee.moves.DESnookerMove(),.2)])
        state=walk;completed=0;history=[];diagnostics=dict(evaluated=False,passed=False);status=None
        while True:
            destination=PILOT_STEPS if completed<PILOT_STEPS else completed+CHECK_EVERY
            if mode=='quick':destination=min(destination,steps)
            if max_steps is not None:destination=min(destination,max_steps)
            made=0
            while completed+made<destination:
                if stop_event is not None and stop_event.is_set():break
                count=min(100,destination-completed-made)
                state=sampler.run_mcmc(state,count,progress=False);made+=count
            if made:
                block=sampler.get_chain();lp=sampler.get_log_prob();store.append(block,lp,state);sampler.reset();completed=store.steps
            if completed==PILOT_STEPS and made:
                pilot,_=store.read();proposal,pinfo=freeze_proposal(pilot,target,best,alt,widths)
                write_json(out/'frozen_proposal.json',pinfo)
                sampler=emcee.EnsembleSampler(walkers,target.ndim,target.log_prob,vectorize=True,
                    moves=[(emcee.moves.DEMove(),.5),(emcee.moves.DESnookerMove(),.2),(emcee.moves.MHMove(proposal),.3)])
                sampler.random_state=state.random_state
            q=None
            if completed>burn:
                chain,lps=store.read();q=(100*np.quantile(chain[burn:,:,2],[.025,.5,.975])).tolist()
                if completed>=burn+2000 and completed%CHECK_EVERY==0:
                    diagnostics=convergence_diagnostics(chain,lps,burn,history);history.append(diagnostics)
                    write_json(out/f'checkpoint_{completed:07d}.json',dict(steps=completed,Ce3_quantiles_percent=q,
                        elapsed_seconds=time.time()-start,convergence=diagnostics))
                del chain,lps
            elapsed=time.time()-start
            msg=f'{completed} steps | elapsed {elapsed/60:.1f} min'
            if q is not None:msg+=f' | Ce³⁺ {q[1]:.3f}% (95% interval {q[0]:.3f}–{q[2]:.3f}%)'
            else:msg+=' | burn-in stage'
            if mode=='convergence' and diagnostics.get('evaluated'):msg+=' | '+('engineering convergence criteria met' if diagnostics['passed'] else 'checking convergence')
            _emit(log,msg)
            event=dict(event='checkpoint',steps=completed,total_steps=steps if mode=='quick' else None,
                elapsed_seconds=elapsed,Ce3_quantiles_percent=q,convergence_checks=diagnostics,
                file_index=file_index,file_count=file_count)
            write_json(out/'progress.json',event);_progress(progress,**event)
            if stop_event is not None and stop_event.is_set():status='cancelled';break
            if mode=='quick' and completed>=steps:status='quick_complete';break
            if mode=='convergence' and diagnostics.get('passed'):status='converged';break
            if max_steps is not None and completed>=max_steps:status='budget_exhausted';break
        result=_export(out,spectrum,target,store,burn,status,config,diagnostics,time.time()-start)
        label={'converged':'Engineering convergence criteria met','quick_complete':'Quick fit complete','cancelled':'Stopped; available results saved','budget_exhausted':'Explicit step budget reached; convergence criteria have not all been met'}[status]
        _emit(log,f'{label}. Output: {out}')
        final_event=dict(event='complete',status=status,steps=store.steps,elapsed_seconds=time.time()-start,
            output_dir=str(out),file_index=file_index,file_count=file_count,
            Ce3_quantiles_percent=[result['Ce3_percent'][k] for k in ('q025','median','q975')] if result.get('posterior_available') else None)
        write_json(out/'progress.json',final_event)
        _progress(progress,**final_event)
        organize_outputs(out)
        return out,status
    except CancelledBeforeSampling:
        stopped=dict(version=__version__,status='cancelled',posterior_available=False,
            total_steps=0,configuration=config,experiment_metadata_file='experiment_metadata.json')
        write_json(out/'fit_result.json',stopped)
        write_result_report(out,stopped)
        organize_outputs(out)
        _emit(log,f'Stopped during initialization. Available files: {out}');return out,'cancelled'
    except Exception as exc:
        exc.run_dir=str(out)
        failure=dict(version=__version__,status='failed',posterior_available=False,error=str(exc),configuration=config,
            total_steps=store.steps if store else 0,experiment_metadata_file='experiment_metadata.json')
        write_json(out/'fit_result.json',failure)
        (out/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        try:
            write_result_report(out,failure)
            organize_outputs(out)
        except Exception:
            # Preserve the original failure even if the secondary report cannot be written.
            pass
        raise

def run_fit(input_paths,*,outdir,model='LA',satellites=False,mode='quick',walkers=80,steps=8000,burn=4000,
            seed=42,background_mode='fixed_shirley',vamas_block=None,vamas_ordinate=None,kinetic_offset_ev=None,
            satellite_gamma=.8,max_steps=None,log=None,progress=None,stop_event=None):
    paths=[Path(p) for p in input_paths]
    if not paths:raise ValueError('Select at least one input file.')
    if model not in ('LA','Voigt') or mode not in ('quick','convergence'):raise ValueError('Unknown model or run mode.')
    if not isinstance(walkers,int) or walkers<22:raise ValueError('At least 22 walkers are required.')
    if burn<PILOT_STEPS or burn%CHECK_EVERY:raise ValueError('Burn-in must be a multiple of 2000 and at least 2000.')
    if mode=='quick' and (steps<=burn or steps%CHECK_EVERY):raise ValueError('Quick-fit steps must exceed burn-in and be a multiple of 2000.')
    if not 0<=seed<2**32:raise ValueError('Random seed must be an integer from 0 through 2**32-1.')
    if max_steps is not None and (max_steps<=burn or max_steps%CHECK_EVERY):raise ValueError('Explicit max_steps must exceed burn-in and be a multiple of 2000.')
    for path in paths:
        if not path.is_file():raise FileNotFoundError(path)
    root=Path(outdir);root.mkdir(parents=True,exist_ok=True)
    if len(paths)>1:root=root/f"batch_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}";root.mkdir()
    results=[]
    for i,path in enumerate(paths,1):
        if len(paths)>1 and stop_event is not None and stop_event.is_set():break
        try:
            output,status=_fit_one(path,root,model=model,satellites=satellites,mode=mode,walkers=walkers,
                steps=steps,burn=burn,seed=seed,background_mode=background_mode,vamas_block=vamas_block,
                vamas_ordinate=vamas_ordinate,kinetic_offset_ev=kinetic_offset_ev,satellite_gamma=satellite_gamma,
                max_steps=max_steps,log=log,progress=progress,stop_event=stop_event,file_index=i,file_count=len(paths))
            results.append(dict(source=str(path),output_dir=str(output),status=status))
        except Exception as exc:
            if len(paths)==1:raise
            results.append(dict(source=str(path),status='failed',error=str(exc),output_dir=getattr(exc,'run_dir',None)));_emit(log,f'{path.name} failed: {exc}')
    if len(paths)==1:return output
    status='cancelled' if stop_event is not None and stop_event.is_set() else 'completed_with_errors' if any(r['status']=='failed' for r in results) else 'batch_complete'
    write_json(root/'batch_summary.json',dict(version=__version__,status=status,requested_files=len(paths),results=results))
    write_batch_report(root,results,status)
    organize_outputs(root)
    return root

def main():
    if sys.stdout is not None and hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    p=argparse.ArgumentParser(description='Ce 3d XPS Fitter 2.0: constrained Bayesian fitting.')
    p.add_argument('--input',nargs='+',required=True);p.add_argument('--outdir',required=True)
    p.add_argument('--model',choices=['LA','Voigt'],default='LA');p.add_argument('--satellites',action='store_true')
    p.add_argument('--mode',choices=['quick','convergence'],default='quick')
    p.add_argument('--walkers',type=int,default=80);p.add_argument('--steps',type=int,default=8000)
    p.add_argument('--burn',type=int,default=4000);p.add_argument('--seed',type=int,default=42)
    p.add_argument('--background-mode',choices=['fixed_shirley','provided'],default='fixed_shirley')
    p.add_argument('--vamas-block',type=int);p.add_argument('--vamas-ordinate',type=int);p.add_argument('--kinetic-offset-ev',type=float)
    p.add_argument('--satellite-gamma',type=float,default=.8)
    p.add_argument('--max-steps',type=int,help='Optional explicit budget for convergence validation; omitted means no fixed cap.')
    p.add_argument('--version',action='version',version=__version__)
    a=p.parse_args();stop=Event()
    try:
        result=run_fit(a.input,outdir=a.outdir,model=a.model,satellites=a.satellites,mode=a.mode,
            walkers=a.walkers,steps=a.steps,burn=a.burn,seed=a.seed,background_mode=a.background_mode,
            vamas_block=a.vamas_block,vamas_ordinate=a.vamas_ordinate,kinetic_offset_ev=a.kinetic_offset_ev,
            satellite_gamma=a.satellite_gamma,max_steps=a.max_steps,log=print,stop_event=stop)
        print(result)
    except KeyboardInterrupt:
        print('Interrupted. Closed chain chunks and metadata remain in the output directory.');raise SystemExit(130)

if __name__=='__main__':main()
