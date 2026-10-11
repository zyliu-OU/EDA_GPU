from dataclasses import replace
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest
from eda_gpu.engine import Settings
from eda_gpu.runtime import WorkloadPacer

class Clock:
    def __init__(self):self.now=0.
    def time(self):return self.now
    def sleep(self,value):self.now+=value


def test_pacing_does_not_double_count_scf_callbacks():
    clock=Clock();p=WorkloadPacer(Settings(compute_duty_cycle=.5),clock=clock.time,sleep=clock.sleep)
    def workload():
        clock.now+=4;p.checkpoint()
        clock.now+=2
        return 7
    assert p.run(workload)==7
    assert p.work_seconds==6 and p.idle_seconds==6 and clock.now==12


def test_temperature_cooldown_hysteresis_and_timeout(monkeypatch):
    clock=Clock();p=WorkloadPacer(Settings(gpu_temperature_limit=75),clock=clock.time,sleep=clock.sleep)
    readings=iter([(76,None),(70,None),(64,None)])
    monkeypatch.setattr(p,'temperatures',lambda:next(readings))
    p.guard();assert p.cooldown_seconds==4
    p=WorkloadPacer(Settings(gpu_temperature_limit=75,max_cooldown_seconds=3),clock=clock.time,sleep=clock.sleep)
    monkeypatch.setattr(p,'temperatures',lambda:(80,None))
    with pytest.raises(RuntimeError,match='cooldown timed out'):p.guard()


def test_invalid_compute_limits():
    with pytest.raises(ValueError,match='duty cycle'):
        WorkloadPacer(Settings(compute_duty_cycle=0))
    with pytest.raises(ValueError,match='resume'):
        WorkloadPacer(Settings(gpu_temperature_limit=60,gpu_resume_temperature=65))


def test_cpu_cap_helper_restores_saved_values_and_preserves_later_changes(tmp_path,monkeypatch):
    path=Path(__file__).resolve().parents[1]/'scripts/cpu_frequency_cap.py'
    spec=importlib.util.spec_from_file_location('cpu_cap_test',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    m.BASE=tmp_path/'policies';m.BASE.mkdir();m.STATE=tmp_path/'state.json';m.LOG=tmp_path/'log.txt'
    for i in range(2):
        policy=m.BASE/f'policy{i}';policy.mkdir()
        for name,value in dict(cpuinfo_min_freq=1000000,cpuinfo_max_freq=5000000,scaling_min_freq=1700000,scaling_max_freq=5000000).items():
            (policy/name).write_text(str(value))
    monkeypatch.setattr(m.subprocess,'Popen',lambda *a,**k:SimpleNamespace(pid=123))
    result=m.apply(3000,60)
    assert result['automatic_restore_minutes']==60
    assert all(m.read(p,'scaling_max_freq')==3000000 for p in m.BASE.iterdir())
    # A later user adjustment must not be overwritten by our automatic restore.
    m.write(m.BASE/'policy1','scaling_max_freq',2500000)
    result=m.restore()
    assert m.read(m.BASE/'policy0','scaling_max_freq')==5000000
    assert m.read(m.BASE/'policy1','scaling_max_freq')==2500000
    assert result['later_changes_preserved']==['policy1'] and not m.STATE.exists()
