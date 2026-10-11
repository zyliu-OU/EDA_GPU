"""Process-scoped workload pacing; never changes global GPU power/clocks."""
import math
import os
import subprocess
import time
from pathlib import Path


def cpu_temperature():
    readings=[]
    for sensor in Path('/sys/class/hwmon').glob('hwmon*'):
        try:
            if sensor.joinpath('name').read_text().strip() not in ('k10temp','coretemp','zenpower'):
                continue
            for path in sensor.glob('temp*_input'):
                value=float(path.read_text())/1000
                if 0<value<130:readings.append(value)
        except (OSError,ValueError):
            continue
    return max(readings) if readings else None


class WorkloadPacer:
    def __init__(self,settings,device=None,clock=time.monotonic,sleep=time.sleep):
        duty=settings.compute_duty_cycle
        if not math.isfinite(duty) or not 0<duty<=1:
            raise ValueError('compute duty cycle must be in (0,1]')
        for high,low in ((settings.gpu_temperature_limit,settings.gpu_resume_temperature),
                         (settings.cpu_temperature_limit,settings.cpu_resume_temperature)):
            if not math.isfinite(high) or not math.isfinite(low) or high<0 or (high and not 0<low<high):
                raise ValueError('temperature limits require 0 < resume < maximum, or maximum=0 to disable')
        if not math.isfinite(settings.max_cooldown_seconds) or settings.max_cooldown_seconds<=0:
            raise ValueError('cooldown timeout must be positive')
        self.settings=settings;self.device=device;self.clock=clock;self.sleep=sleep
        self.work_seconds=0.;self.idle_seconds=0.;self.cooldown_seconds=0.;self.last=clock()
        self.max_gpu_temperature=None;self.max_cpu_temperature=None

    def temperatures(self):
        gpu=None
        if self.device is not None and self.settings.gpu_temperature_limit:
            result=subprocess.run(['nvidia-smi','-i',self.device,'--query-gpu=temperature.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=10)
            try:
                gpu=float(result.stdout.strip())
                if result.returncode or not math.isfinite(gpu):raise ValueError('invalid reading')
            except ValueError as exc:
                raise RuntimeError('GPU temperature guard cannot read the selected GPU sensor') from exc
            self.max_gpu_temperature=max(gpu,self.max_gpu_temperature or gpu)
        cpu=cpu_temperature() if self.settings.cpu_temperature_limit else None
        if cpu is not None:self.max_cpu_temperature=max(cpu,self.max_cpu_temperature or cpu)
        if getattr(self,"progress",None):
            self.progress.update(self.progress.current.get("stage","temperature check"),gpu_temperature_c=gpu,cpu_temperature_c=cpu)
        return gpu,cpu

    def pause(self,seconds,cooldown=False):
        if seconds<=0:return
        begin=self.clock();self.sleep(seconds);elapsed=self.clock()-begin
        self.idle_seconds+=elapsed
        if cooldown:self.cooldown_seconds+=elapsed

    def guard(self):
        gpu,cpu=self.temperatures()
        hot_gpu=gpu is not None and gpu>=self.settings.gpu_temperature_limit
        hot_cpu=cpu is not None and cpu>=self.settings.cpu_temperature_limit
        if not (hot_gpu or hot_cpu):return
        if getattr(self,"progress",None):self.progress.update("temperature cooldown",gpu_temperature_c=gpu,cpu_temperature_c=cpu)
        deadline=self.clock()+self.settings.max_cooldown_seconds
        while ((hot_gpu and (gpu is None or gpu>self.settings.gpu_resume_temperature)) or
               (hot_cpu and (cpu is None or cpu>self.settings.cpu_resume_temperature))):
            if self.clock()>=deadline:raise RuntimeError('temperature cooldown timed out; calculation stopped')
            self.pause(min(2.,deadline-self.clock()),cooldown=True)
            gpu,cpu=self.temperatures()

    def pace(self,work):
        work=max(0.,work);self.work_seconds+=work
        self.pause(work*(1/self.settings.compute_duty_cycle-1))
        self.guard();self.last=self.clock()

    def checkpoint(self,env=None):
        self.pace(self.clock()-self.last)

    def run(self,function,*args,**kwargs):
        self.guard();self.last=self.clock()
        begin=self.last;work=self.work_seconds;idle=self.idle_seconds
        result=function(*args,**kwargs)
        # SCF callbacks have already paced their work. Only pace the remainder.
        self.pace(self.clock()-begin-(self.idle_seconds-idle)-(self.work_seconds-work))
        return result

    def report(self):
        return dict(scope='this calculation process only; GPU power and clocks unchanged',
            cpu_thread_limit=min(self.settings.threads,self.settings.cpu_core_limit) if self.settings.cpu_core_limit else self.settings.threads,
            allowed_cpu_cores=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None,
            duty_cycle_target=self.settings.compute_duty_cycle,work_seconds=self.work_seconds,
            pacing_and_cooldown_seconds=self.idle_seconds,temperature_cooldown_seconds=self.cooldown_seconds,
            max_gpu_temperature_c=self.max_gpu_temperature,max_cpu_temperature_c=self.max_cpu_temperature,
            cpu_sensor_available=self.max_cpu_temperature is not None,
            gpu_temperature_limit_c=self.settings.gpu_temperature_limit,
            cpu_temperature_limit_c=self.settings.cpu_temperature_limit,
            limitation='checked at SCF-cycle and workload boundaries; long kernels are not interrupted; no per-process watt or instantaneous-utilization guarantee')
