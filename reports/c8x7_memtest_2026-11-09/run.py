import json, sys, time
from pathlib import Path
import numpy as np
from eda_gpu.engine import Settings, singlepoint
from eda_gpu.model import load
S = Path(sys.argv[1])
symbols, coords, meta = load(S/'seed/patch_tributylammonium_0.xyz', S/'seed/patch_tributylammonium_0.json')
st = Settings(backend='gpu', xc='m06-2x', basis='def2-svp', max_cycle=200, threads=16, solvent='cpcm',
              solvent_eps=78.3553, compute_duty_cycle=1.0, gpu_temperature_limit=80.0, cpu_temperature_limit=85.0)
r = singlepoint(symbols, coords, meta, st, out=S/'progress')
r.pop('gradient_hartree_per_bohr_full', None)
g = np.array(r['gradient_hartree_per_bohr'])
print(json.dumps(dict(energy=r['energy_hartree'], gmax=float(abs(g).max()), seconds=r['seconds'], limits=r['compute_limits']), indent=1))
