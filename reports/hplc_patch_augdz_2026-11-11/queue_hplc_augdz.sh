#!/bin/bash
# aug-cc-pVDZ GKS-EDA(sol) single points on the lowest-energy patch placement per reagent (def2-SVP geometries).
# Direct integrals (no DF tensor); smallest case first as the memory/timing test. Resumable.
cd /home/zyliu/Documents/EDA_GPU
source scripts/gpu-env.sh
O=runs/hplc_patch_augdz; mkdir -p $O
( while true; do echo "$(date +%T) $(nvidia-smi --query-gpu=memory.used,temperature.gpu --format=csv,noheader,nounits) $(free -m | awk '/Mem:/{print $3}') $(free -m | awk '/Swap:/{print $3}')"; sleep 10; done ) > $O/mem.log &
MON=$!
R="runs/hplc_patch_cpcm_def2svp runs/hplc_patch_sampling_cpcm_def2svp"
for c in patch_teah_2 patch_hexylammonium_0 patch_dipeah_0 patch_dibutylammonium_0 patch_tributylammonium_3; do
  .env/bin/python scripts/sensitivity_hplc.py $R --variants augdz --cases $c --output $O
done
kill $MON
echo "$(date -Is) queue done"
