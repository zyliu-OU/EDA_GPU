while true; do
  echo "$(date +%T) $(nvidia-smi --query-gpu=memory.used,temperature.gpu --format=csv,noheader,nounits) $(free -m | awk '/Mem:/{print $3}') $(free -m | awk '/Swap:/{print $3}')"
  sleep 5
done
