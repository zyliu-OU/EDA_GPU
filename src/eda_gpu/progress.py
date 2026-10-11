"""Atomic run snapshots and read-only terminal monitoring."""
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import time


def atomic_json(path,data):
    path=Path(path);temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temporary.replace(path)


class ProgressWriter:
    def __init__(self,out):
        self.out=Path(out);self.started=time.monotonic();self.sequence=0
        from .citations import information
        self.current=dict(**information(),schema_version=1,pid=os.getpid(),status='running')
    def update(self,stage,**values):
        self.sequence+=1
        event=dict(timestamp=datetime.now(timezone.utc).isoformat(),sequence=self.sequence,
                   elapsed_seconds=time.monotonic()-self.started,stage=stage,**values)
        self.current.update(event)
        atomic_json(self.out/'progress.json',self.current)
        with (self.out/'events.jsonl').open('a') as stream:
            stream.write(json.dumps(event,allow_nan=False)+'\n')
    def scf(self,env,label='SCF'):
        values=dict(scf_cycle=int(env.get('cycle',-1))+1)
        for source,target in [('e_tot','energy_hartree'),('norm_gorb','orbital_gradient_norm'),('norm_ddm','density_change_norm')]:
            value=env.get(source)
            if value is not None:
                value=float(value.item() if hasattr(value,'item') else value)
                if math.isfinite(value):values[target]=value
        self.update(label,**values)


def snapshot(out,log=None):
    from .citations import information
    out=Path(out);result=dict(**information(),run_folder=str(out),status='waiting',stage='awaiting progress')
    for filename in ('progress.json','result.json'):
        try:
            data=json.loads((out/filename).read_text())
        except (OSError,ValueError):continue
        if filename=='progress.json':result.update(data)
        else:
            result.update(information(data.get('metadata')))
            for key in ('research_acknowledgement','project_reference','method_references'):
                if key in data:result[key]=data[key]
            result.update(status=data.get('status','complete'),stage='finished',result_file=str(out/filename))
            for key in ('error','energy_hartree','optimization_converged','components_hartree','counterpoise_interaction_energy_hartree','total_closure_residual_hartree','backend_used','scientifically_validated'):
                if key in data:result[key]=data[key]
    if log:
        try:
            with Path(log).open('rb') as stream:
                stream.seek(0,2);length=stream.tell();stream.seek(max(0,length-65536));lines=stream.read().decode(errors='replace').splitlines()
            lines=[re.sub(r'\x1b\[[0-9;]*m','',line) for line in lines]
            steps=[line for line in lines if line.startswith('Step ')]
            result['latest_log_lines']=steps[-2:] if steps else lines[-4:]
            if result['status']=='waiting' and lines:result.update(status='running',stage='following legacy calculation log')
        except OSError:result['log_status']='not available yet'
    return result


def render(data,include_citations=True):
    lines=[f"EDA_GPU | {data['status']} | {data['stage']}",f"Run: {data['run_folder']}"]
    if 'scf_cycle' in data:lines.append(f"SCF cycle: {data['scf_cycle']}")
    if 'optimization_evaluation' in data:lines.append(f"Optimization evaluation: {data['optimization_evaluation']}")
    if 'energy_hartree' in data:lines.append(f"Energy: {data['energy_hartree']:.12f} hartree")
    for key in ('gradient_rms_hartree_per_bohr','gradient_max_hartree_per_bohr'):
        if key in data:lines.append(f"{key}: {data[key]:.6e}")
    for key in ('gpu_temperature_c','cpu_temperature_c'):
        if data.get(key) is not None:lines.append(f"{key}: {data[key]:.1f} °C")
    if 'components_hartree' in data:
        lines.append('GKS-EDA components (hartree):')
        lines.extend(f"  {key}: {value:+.12f}" for key,value in data['components_hartree'].items())
    if 'counterpoise_interaction_energy_hartree' in data:
        lines.append(f"CP interaction: {data['counterpoise_interaction_energy_hartree']:+.12f} hartree")
    if data.get('scientifically_validated') is False:lines.append('Experimental decomposition; external scientific validation pending.')
    if 'error' in data:lines.append('Error: '+data['error'])
    lines.extend(data.get('latest_log_lines',[]))
    from .citations import text
    if include_citations:lines.append(text(data))
    return '\n'.join(lines)


def monitor(out,log=None,interval=1.,once=False,json_output=False):
    if not math.isfinite(interval) or interval<=0:raise ValueError('monitor interval must be positive')
    previous=None
    try:
        while True:
            data=snapshot(out,log);serialized=json.dumps(data,sort_keys=True)
            if serialized!=previous:
                print(json.dumps(data,indent=2) if json_output else render(data,include_citations=previous is None),flush=True)
                previous=serialized
            if once or data['status'] in ('success','complete','failed','not_converged'):
                return 1 if data['status'] in ('failed','not_converged') else 0
            time.sleep(interval)
    except KeyboardInterrupt:
        return 0
