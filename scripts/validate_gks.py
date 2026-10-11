#!/usr/bin/env python
"""Component-wise CPU/GPU checks and repeated complete experimental EDA timing."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import resource
import statistics
import time
import numpy as np
from eda_gpu.engine import Settings, hardware
from eda_gpu.gks import Algebra, decompose
from eda_gpu.model import load

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--gpu',action='store_true',help='require actual GPU EDA; never substitutes CPU')
    parser.add_argument('--repeats',type=int,default=3)
    args=parser.parse_args()
    if args.repeats<1: parser.error('--repeats must be positive')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    result={'hardware':hardware(),'gpu_requested':args.gpu,'gpu_executed':False,
            'scientifically_validated':False,'external_reference_validation':'not established',
            'component_agreement_tolerance_hartree':1e-7,'closure_tolerance_hartree':1e-8,
            'warmup_policy':'one full decomposition per system/backend excluded; fresh reference SCFs on every repetition',
            'memory_note':'process cumulative RSS and CuPy pool/device snapshots, not peak EDA GPU memory',
            'cases':[],'passed':False}
    try:
        backends=['cpu','gpu'] if args.gpu else ['cpu']
        cases=[('water_dimer','def2-svp',False),('co2_dimer_xeda_ii','6-31g*',True)]
        for name,basis,cart in cases:
            symbols,coords,meta=load(ROOT/'examples'/f'{name}.xyz',ROOT/'examples'/f'{name}.json')
            case={'system':name,'basis':basis,'cart':cart,'backends':{}}
            first={}
            for backend in backends:
                algebra=Algebra(backend)
                settings=Settings(backend=backend,basis=basis,cart=cart,density_fitting=False,conv_tol=1e-12)
                durations=[];measurements=[]
                for repetition in range(args.repeats+1):
                    run=out/f'{name}_{backend}_{repetition}';run.mkdir()
                    algebra.synchronize();start=time.perf_counter()
                    calc=decompose(symbols,coords,meta,settings,run)
                    algebra.synchronize();elapsed=time.perf_counter()-start
                    (run/'result.json').write_text(json.dumps(calc,indent=2,allow_nan=False)+'\n')
                    if backend=='gpu': result['gpu_executed']=True
                    if repetition==0: first[backend]=calc
                    else:
                        durations.append(elapsed)
                        measurements.append({'seconds':elapsed,'calculation_seconds':calc['calculation_seconds'],
                            'components_hartree':calc['components_hartree'],'closure_residual_hartree':calc['total_closure_residual_hartree']})
                    if abs(calc['total_closure_residual_hartree'])>result['closure_tolerance_hartree']:
                        raise RuntimeError('closure exceeds threshold')
                    print(f'{name} {backend} repetition {repetition}: {elapsed:.3f}s',flush=True)
                data={'settings':asdict(settings),'repetitions':measurements,'median_seconds':statistics.median(durations),
                      'process_cumulative_max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
                if backend=='gpu':
                    xp=algebra.xp; free,total=xp.cuda.runtime.memGetInfo()
                    data['gpu_memory_snapshot_bytes']={'pool_used':xp.get_default_memory_pool().used_bytes(),
                        'pool_total':xp.get_default_memory_pool().total_bytes(),'device_free':free,'device_total':total}
                case['backends'][backend]=data
            if args.gpu:
                # Compare every recorded GPU repetition to the CPU baseline; no total-only cancellation.
                errors={k:max(abs(m['components_hartree'][k]-v) for m in case['backends']['gpu']['repetitions'])
                        for k,v in first['cpu']['components_hartree'].items()}
                case['cpu_gpu_max_component_difference_hartree']=errors
                case['cpu_gpu_component_agreement_passed']=max(errors.values())<result['component_agreement_tolerance_hartree']
                if not case['cpu_gpu_component_agreement_passed']: raise RuntimeError(f'{name}: component mismatch {errors}')
            if name=='co2_dimer_xeda_ii':
                reference=json.loads((ROOT/'references/xeda_ii_co2_b3lyp.json').read_text())
                values=first['cpu']['components_hartree']
                diffs={k:values[k]-v for k,v in reference['components'].items() if k in values}
                diffs['total']=first['cpu']['counterpoise_interaction_energy_hartree']-reference['components']['total']
                case['published_reference']={'source_url':reference['source_url'],'computed_minus_published_hartree':diffs,
                    'status':'exploratory; incomplete specification of external settings','reported_precision_hartree':reference['reported_precision_hartree'],
                    'unknown_settings':reference['unspecified_parameters']}
            result['cases'].append(case)
        result['passed']=True
    except Exception as exc:
        result['error']=f'{type(exc).__name__}: {exc}'
    (out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':result['passed'],'result':str(out/'result.json'),'error':result.get('error')},indent=2))
    return 0 if result['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
