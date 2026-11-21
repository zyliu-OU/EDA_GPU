#!/usr/bin/env bash
cd /home/zyliu/Documents/EDA_GPU; source scripts/gpu-env.sh; D=reports/c8x7_memtest_2026-11-09
bash $D/mon.sh > $D/mem.log 2>&1 & MON=$!
/usr/bin/time -v .env/bin/python $D/run.py $D > $D/run.out 2> $D/run.err; echo "rc=$?" >> $D/run.out
kill $MON
