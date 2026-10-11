#!/usr/bin/python3
"""Temporary shared CPU frequency cap. Requires root; never persists over reboot.
Usage: sudo /usr/bin/python3 scripts/cpu_frequency_cap.py apply --max-mhz 3000 --minutes 60
       sudo /usr/bin/python3 scripts/cpu_frequency_cap.py restore
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid

STATE=Path('/run/eda_gpu_cpu_frequency_cap.json')
LOG=Path('/run/eda_gpu_cpu_frequency_restore.log')
BASE=Path('/sys/devices/system/cpu/cpufreq')

def read(policy,key):return int((policy/key).read_text())
def write(policy,key,value):(policy/key).write_text(str(value))

def restore(token=None):
    if not STATE.exists():return {'status':'already_restored'}
    state=json.loads(STATE.read_text())
    if token and token!=state['token']:return {'status':'different_cap_active'}
    restored=[];skipped=[]
    for item in state['policies']:
        name=item['name']
        if not re.fullmatch(r'policy[0-9]+',name):raise RuntimeError('invalid saved policy')
        policy=BASE/name
        current=read(policy,'scaling_max_freq')
        # Respect a later manual/system change instead of overwriting it.
        if current not in (item['applied_max'],item['original_max']):
            skipped.append(name);continue
        write(policy,'scaling_max_freq',item['original_max'])
        if read(policy,'scaling_min_freq')==item['applied_min']:
            write(policy,'scaling_min_freq',item['original_min'])
        restored.append(name)
    STATE.unlink()
    return {'status':'restored','policies':restored,'later_changes_preserved':skipped}

def apply(mhz,minutes):
    if not 0<mhz<10000 or not 0<minutes<=240:raise ValueError('invalid frequency or duration')
    if STATE.exists():raise RuntimeError('a saved cap already exists; restore it before applying another')
    policies=sorted(BASE.glob('policy*'))
    if not policies:raise RuntimeError('CPU frequency control unavailable')
    target=mhz*1000;items=[]
    for policy in policies:
        if not read(policy,'cpuinfo_min_freq')<=target<=read(policy,'cpuinfo_max_freq'):
            raise ValueError('requested frequency is outside supported policy range')
        original_min=read(policy,'scaling_min_freq');original_max=read(policy,'scaling_max_freq')
        items.append(dict(name=policy.name,original_min=original_min,original_max=original_max,
                          applied_min=min(original_min,target),applied_max=min(original_max,target)))
    token=uuid.uuid4().hex
    state=dict(token=token,scope='shared CPU hardware; all applications affected temporarily',policies=items,restore_after_seconds=minutes*60)
    fd=os.open(STATE,os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW,0o600)
    with os.fdopen(fd,'w') as stream:json.dump(state,stream,indent=2)
    try:
        for item in items:
            policy=BASE/item['name'];write(policy,'scaling_min_freq',item['applied_min']);write(policy,'scaling_max_freq',item['applied_max'])
            item['applied_max']=read(policy,'scaling_max_freq')
        STATE.write_text(json.dumps(state,indent=2))
        log_fd=os.open(LOG,os.O_WRONLY|os.O_CREAT|os.O_APPEND|os.O_NOFOLLOW,0o600)
        with os.fdopen(log_fd,'a') as log:
            child=subprocess.Popen(['/usr/bin/python3',str(Path(__file__).resolve()),'_timer','--token',token,'--seconds',str(minutes*60)],stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True,close_fds=True)
    except Exception:
        restore();raise
    return dict(status='applied',max_mhz=mhz,policy_count=len(items),automatic_restore_minutes=minutes,restore_process_pid=child.pid,scope=state['scope'])

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('apply');a.add_argument('--max-mhz',type=int,default=3000);a.add_argument('--minutes',type=int,default=60)
    sub.add_parser('restore');sub.add_parser('status')
    t=sub.add_parser('_timer');t.add_argument('--token',required=True);t.add_argument('--seconds',type=int,required=True)
    args=p.parse_args()
    if args.command=='status':print(STATE.read_text() if STATE.exists() else '{"status":"no saved cap"}');return
    if os.geteuid()!=0:raise SystemExit('Administrator access required: run this helper with sudo in your terminal.')
    if args.command=='_timer':time.sleep(args.seconds);result=restore(args.token)
    elif args.command=='restore':result=restore()
    else:result=apply(args.max_mhz,args.minutes)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
