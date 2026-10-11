#!/usr/bin/env bash
# Source before using the isolated CUDA 12 wheel installation.
eda_gpu_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
eda_gpu_site="$eda_gpu_root/.env/lib/python3.12/site-packages"
for eda_gpu_lib in "$eda_gpu_site"/nvidia/*/lib "$eda_gpu_site"/cutensor/lib; do
    if [[ -d "$eda_gpu_lib" ]]; then
        export LD_LIBRARY_PATH="$eda_gpu_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    fi
done
unset eda_gpu_root eda_gpu_site eda_gpu_lib
