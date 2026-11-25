"""Adaptive piecewise Shirley anchored at the lowest measured crossing point.

    This is preprocessing, not a jointly inferred background. No pointwise
    clipping is used. Exact non-crossing can insert many noise-driven anchors.
"""
import numpy as np
from smart_like_shirley_background import shirley_standard

def shirley_v7(energy, intensity, *, initial_anchor_indices=None,
               edge_fraction=.10, crossing_tolerance=0., max_refinements=None,
               smoothing_window=1, endpoint_mode='edge_min'):
    e=np.asarray(energy,float); y=np.asarray(intensity,float)
    if e.ndim!=1 or y.shape!=e.shape or len(e)<12:
        raise ValueError('Require equal 1D arrays with at least 12 points')
    if not np.all(np.isfinite(e)) or not np.all(np.isfinite(y)):
        raise ValueError('Nonfinite data')
    if crossing_tolerance<0:raise ValueError('Negative crossing tolerance')
    order=np.argsort(e);e=e[order];y=y[order]
    if np.any(np.diff(e)<=0):raise ValueError('Duplicate energies')
    raw=y.copy()
    if not isinstance(smoothing_window,int) or smoothing_window<1 or smoothing_window%2!=1:
        raise ValueError('Smoothing window must be a positive odd integer')
    if smoothing_window>len(y):raise ValueError('Smoothing window exceeds data length')
    if smoothing_window>1:
        half=smoothing_window//2
        # Truncated local means: no padding or invented exterior intensities.
        y=np.array([raw[max(0,i-half):min(len(raw),i+half+1)].mean() for i in range(len(raw))])
    n=len(e); m=max(6,int(n*edge_fraction))
    if not 0<edge_fraction<=.5:raise ValueError('Invalid edge fraction')
    if initial_anchor_indices is None:
        if endpoint_mode=='ends':left,right=0,n-1
        elif endpoint_mode=='edge_min':
            left=int(np.argmin(y[:m]));right=n-m+int(np.argmin(y[-m:]))
        else:raise ValueError('Unknown endpoint mode')
        anchors={left,right}
    else:
        supplied=np.asarray(initial_anchor_indices)
        if not np.issubdtype(supplied.dtype,np.integer) or np.any((supplied<0)|(supplied>=n)):
            raise ValueError('Anchor indices must address original input rows')
        inverse=np.empty(n,int);inverse[order]=np.arange(n)
        anchors=set(inverse[supplied].tolist())
        if len(anchors)<2:raise ValueError('At least two anchors required')
    # Constant exterior baseline, then refine it too if a crossing is found.
    limit=n if max_refinements is None else int(max_refinements)
    if limit<0:raise ValueError('Negative refinement limit')
    history=[]
    def construct():
        a=sorted(anchors);b=np.empty(n)
        b[:a[0]+1]=y[a[0]];b[a[-1]:]=y[a[-1]]
        for l,r in zip(a[:-1],a[1:]):
            b[l:r+1]=shirley_standard(e[l:r+1],y[l:r+1],y[l],y[r])
        return b
    numerical_eps=32*np.finfo(float).eps*max(1.,float(np.max(abs(y))))
    for iteration in range(limit+1):
        b=construct();cross=np.flatnonzero(b-y>crossing_tolerance+numerical_eps)
        candidates=np.array([i for i in cross if i not in anchors],int)
        if not len(cross) or iteration==limit or not len(candidates):break
        idx=int(candidates[np.argmin(y[candidates])])
        history.append({'iteration':iteration+1,'energy_eV':float(e[idx]),
                        'intensity':float(y[idx]),'original_row_index':int(order[idx]),
                        'crossing_count_before':int(len(cross)),
                        'maximum_crossing_before':float(np.max(b-y))})
        anchors.add(idx)
    original=np.empty(n);original[order]=b
    return original,{'version':'v7_adaptive_piecewise_shirley',
        'converged_no_crossing':bool(not len(cross)),
        'crossing_tolerance':float(crossing_tolerance),'numerical_epsilon':numerical_eps,
        'remaining_crossing_count':int(len(cross)),
        'smoothing_window_points':smoothing_window,'endpoint_mode':endpoint_mode,
        'crossing_detection_target':'local-mean spectrum' if smoothing_window>1 else 'raw spectrum',
        'raw_crossing_count':int(np.sum(b-raw>crossing_tolerance+numerical_eps)),
        'raw_maximum_crossing':float(max(0,np.max(b-raw))),
        'maximum_crossing':float(max(0,np.max(b-y))),
        'anchors':[{'energy_eV':float(e[i]),'intensity':float(y[i]),'original_row_index':int(order[i])} for i in sorted(anchors)],
        'refinements':history,'method':'piecewise re-computation, no pointwise clipping',
        'external_behavior':'constant outside outermost anchors, refined if crossing',
        'warning':'Non-crossing does not establish physical correctness; noise or broad peak valleys can introduce excessive anchors.'}
