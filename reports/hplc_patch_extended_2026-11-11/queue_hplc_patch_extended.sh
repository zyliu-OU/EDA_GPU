#!/bin/bash
# Step 3: BA+, DEA+, TPA+ on the C8 patch, same protocol as the first five reagents (10-27 + 10-28):
# placements 0/1 -> xTB pre-opt + DFT opt + GKS-EDA(sol); xTB screen of 2-9 (--known = 0/1), best 3 -> DFT.
# Waits for the aug-cc-pVDZ queue to free the GPU. Resumable.
cd /home/zyliu/Documents/EDA_GPU
source scripts/gpu-env.sh
while systemctl --user is-active --quiet hplc-augdz; do sleep 300; done
echo "$(date -Is) start"
O=runs/hplc_patch_extended_cpcm_def2svp; S=examples/hplc/patch_extended
G="diethylammonium butylammonium tripropylammonium"
.env/bin/python scripts/run_hplc_stage.py $(for g in $G; do echo $S/patch_${g}_0 $S/patch_${g}_1; done) --preopt xtb --maxsteps 300 --output $O
.env/bin/python scripts/screen_hplc_patch.py $(for g in $G; do for i in 2 3 4 5 6 7 8 9; do echo $S/patch_${g}_$i; done; done) --known $O --pick 3 --output $O/screen
seeds=$(.env/bin/python -c "
import json; s=json.load(open('$O/screen/selection.json'))
print(' '.join(v['picked_seeds'][k] for k in range(3) for v in s.values() if len(v['picked_seeds']) > k))")
.env/bin/python scripts/run_hplc_stage.py $seeds --preopt xtb --maxsteps 300 --output $O
.env/bin/python scripts/nonpolar_hplc.py $O --output $O/nonpolar.json
.env/bin/python scripts/analyze_hplc.py $O --output $O/analysis.json
echo "$(date -Is) queue done"
