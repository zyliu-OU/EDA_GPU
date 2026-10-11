import json
from pathlib import Path
from eda_gpu.progress import ProgressWriter,snapshot,monitor
from eda_gpu.cli import main


def test_live_scf_snapshot_and_completed_components(tmp_path,capsys):
    writer=ProgressWriter(tmp_path)
    writer.scf({'cycle':2,'e_tot':-100.5,'norm_gorb':.001},'fragment 0 SCF')
    data=snapshot(tmp_path)
    assert data['scf_cycle']==3 and data['status']=='running'
    assert data['energy_hartree']==-100.5
    assert len((tmp_path/'events.jsonl').read_text().splitlines())==1
    result=dict(status='success',components_hartree={'electrostatic':-.01},scientifically_validated=False)
    (tmp_path/'result.json').write_text(json.dumps(result))
    assert monitor(tmp_path,interval=.01)==0
    output=capsys.readouterr().out
    assert 'electrostatic' in output and '10.1021/jp500405s' in output


def test_status_legacy_log_and_incomplete_json(tmp_path,capsys):
    (tmp_path/'progress.json').write_text('{unfinished')
    log=tmp_path/'calc.log';log.write_text('Step 5 : Energy = -1.234\n')
    assert main(['status',str(tmp_path),'--log',str(log)])==0
    assert 'Step 5' in capsys.readouterr().out
    assert not (tmp_path/'events.jsonl').exists()
    assert snapshot(tmp_path)['status']=='waiting'


def test_failed_result_status_returns_failure(tmp_path,capsys):
    (tmp_path/'result.json').write_text(json.dumps(dict(status='failed',error='SCF did not converge')))
    assert main(['status',str(tmp_path)])==1
    assert 'SCF did not converge' in capsys.readouterr().out
