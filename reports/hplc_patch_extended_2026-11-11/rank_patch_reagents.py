import json, sys, numpy as np
from scipy.stats import spearmanr
rows=json.load(open(sys.argv[1]))
scr={}
for p in ['runs/hplc_patch_sampling_cpcm_def2svp/screen/screen.json','runs/hplc_patch_extended_cpcm_def2svp/screen/screen.json']:
    scr.update(json.load(open(p))['cases'])
pub=['diethylammonium','teah','dipeah','butylammonium','tripropylammonium','dibutylammonium','hexylammonium','tributylammonium']
kT=0.0019872041*298.15; H=627.509474
out={'published_order_least_to_most_retained':pub,'source':'Donegan Table 4, 100 mM acetate, 20T %ACN; TBA placed last','kT_kcal_mol':kT,'reagents':{},'spearman':{}}
def grp(g): return {k:v for k,v in rows.items() if k.startswith(f'patch_{g}_') and k[len(f'patch_{g}_'):].isdigit()}
for g in pub:
    c=grp(g); E0=min(v['E'] for v in c.values())
    for v in c.values(): v['dE']=(v['E']-E0)*H
    w=np.array([np.exp(-v['dE']/kT) for v in c.values()]); w/=w.sum()
    low=min(c,key=lambda k:c[k]['dE'])
    d=dict(placements=sorted(c), lowest_energy=low, weights=dict(zip(c,w.round(4).tolist())))
    for q in ['tot','epc','np','elec','pauli','pol','gc','des']:
        a=np.array([c[k][q] for k in c])
        d[q]=dict(lowest_energy=c[low][q], best=(a.max() if q=='np' else a.min()), boltzmann=float(w@a), mean=float(a.mean()))
    out['reagents'][g]=d
rk=np.arange(len(pub))
for q in ['tot','epc','np','gc','elec']:
    for m in ['lowest_energy','best','boltzmann','mean']:
        x=[out['reagents'][g][q][m] for g in pub]; x=[-v for v in x] if q!='np' else x
        r,p=spearmanr(rk,x); out['spearman'][f'{q}/{m}']=dict(rho=round(r,2),p=round(p,3),order=[pub[i] for i in np.argsort(x)])
for sub,lab in [(pub[:7],'n7_noTBA'),([g for g in pub if g not in('diethylammonium','butylammonium','tripropylammonium')],'n5_original')]:
    for q in ['tot','epc']:
        for m in ['best','boltzmann']:
            x=[-out['reagents'][g][q][m] for g in sub]; r,p=spearmanr(range(len(sub)),x)
            out['spearman'][f'{lab}:{q}/{m}']=dict(rho=round(r,2),p=round(p,3))
json.dump(rows,open(sys.argv[1],'w'),indent=1); json.dump(out,open(sys.argv[2],'w'),indent=1)
print(f"{'reagent':18s} {'n':>2} {'tot_best':>8} {'tot_Bz':>7} {'tot_mean':>8} {'+CDS_Bz':>7} {'npSASA_Bz':>9} {'gc_Bz':>6} {'elec_Bz':>7}")
for g in pub:
    d=out['reagents'][g]; print(f"{g:18s} {len(d['placements']):2d} {d['tot']['best']:8.2f} {d['tot']['boltzmann']:7.2f} {d['tot']['mean']:8.2f} {d['epc']['boltzmann']:7.2f} {d['np']['boltzmann']:9.0f} {d['gc']['boltzmann']:6.2f} {d['elec']['boltzmann']:7.2f}")
for k,v in out['spearman'].items(): print(k,v['rho'],v['p'],[o[:5] for o in v.get('order',[])])
