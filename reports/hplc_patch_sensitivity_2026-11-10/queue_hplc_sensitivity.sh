#!/bin/bash
# Solvent sensitivity on the 25 C8-patch placements: EDA single points at the water/vdw 1.2 geometries.
# 1) reproduce one reference EDA; 2) vdw_scale 1.4; 3) acetonitrile eps 35.688. Resumable.
cd /home/zyliu/Documents/EDA_GPU
source scripts/gpu-env.sh
O=runs/hplc_patch_sensitivity
R="runs/hplc_patch_cpcm_def2svp runs/hplc_patch_sampling_cpcm_def2svp"
.env/bin/python scripts/sensitivity_hplc.py $R --variants ref --cases patch_teah_2 --output $O
.env/bin/python scripts/sensitivity_hplc.py $R --variants vdw14 acn --output $O
echo "$(date -Is) queue done"
