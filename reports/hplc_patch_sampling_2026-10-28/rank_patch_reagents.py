import json, sys, numpy as np
from scipy.stats import spearmanr
rows=json.load(open(sys.argv[1])); scr=json.load(open('runs/hplc_patch_sampling_cpcm_def2svp/screen/screen.json'))['cases']
pub=['teah','dipeah','dibutylammonium','hexylammonium','tributylammonium']  # least -> most retained
kT=0.0019872041*298.15; H=627.509474
out={'published_order_least_to_most_retained':pub,'kT_kcal_mol':kT,'reagents':{},'spearman':{}}
for g in pub:
    c={k:v for k,v in rows.items() if k.startswith(f'patch_{g}_')}
    E0=min(v['E'] for v in c.values())
    for v in c.values(): v['dE']=(v['E']-E0)*H
    w=np.array([np.exp(-v['dE']/kT) for v in c.values()]); w/=w.sum()
    low=min(c,key=lambda k:c[k]['dE'])
    d=dict(placements=sorted(c), lowest_energy=low, weights=dict(zip(c,w.round(4).tolist())))
    for q in ['tot','epc','np']:
        a=np.array([c[k][q] for k in c])
        d[q]=dict(lowest_energy=c[low][q], best=(a.max() if q=='np' else a.min()), boltzmann=float(w@a), mean=float(a.mean()),
                  first_two=float(min(c[f'patch_{g}_0'][q],c[f'patch_{g}_1'][q])) if q!='np' else float(max(c[f'patch_{g}_0'][q],c[f'patch_{g}_1'][q])))
    out['reagents'][g]=d
rk=np.arange(5)
for q in ['tot','epc','np']:
    for m in ['lowest_energy','best','boltzmann','mean','first_two']:
        x=[out['reagents'][g][q][m] for g in pub]; x=[-v for v in x] if q!='np' else x
        r,p=spearmanr(rk,x); order=[pub[i] for i in np.argsort(x)]
        out['spearman'][f'{q}/{m}']=dict(rho=round(r,2),p=round(p,2),order=order)
# xTB vs DFT on the 15 screened placements
xs=[k for k in rows if k in scr]
xt=[scr[k]['dE_xtb_kcal_mol'] for k in xs]
out['xtb_vs_dft']=dict(n=len(xs), rho_eda=round(spearmanr(xt,[rows[k]['tot'] for k in xs])[0],2), rho_eda_cds=round(spearmanr(xt,[rows[k]['epc'] for k in xs])[0],2),
   within_reagent_top_pick_is_best_dft={g:(min((k for k in xs if k.startswith(f'patch_{g}_')),key=lambda k:scr[k]['dE_xtb_kcal_mol'])==min((k for k in xs if k.startswith(f'patch_{g}_')),key=lambda k:rows[k]['tot'])) for g in pub},
   pairs={k:dict(xtb=round(scr[k]['dE_xtb_kcal_mol'],2),eda=round(rows[k]['tot'],2)) for k in xs})
json.dump(rows,open(sys.argv[1],'w'),indent=1); json.dump(out,open(sys.argv[2],'w'),indent=1)
for g in pub:
    d=out['reagents'][g]; print(g, 'low',d['lowest_energy'], {q:{m:round(v,2) for m,v in d[q].items()} for q in ['tot','epc','np']}, 'w',{k[-1]:x for k,x in d['weights'].items() if x>0.01})
for k,v in out['spearman'].items(): print(k,v)
print({k:v for k,v in out['xtb_vs_dft'].items() if k!='pairs'})
