#!/bin/bash
# Machine state every 15 s: 1-minute load, frequency of the pinned core (CORE, default 2), CPU package
# temperature, available memory, used swap, read rate of the input disk (DISK, default sda).
#   CORE=2 DISK=sda monitor.sh monitor.log
CORE=${CORE:-2}
DISK=${DISK:-sda}
OUT=${1:-monitor.log}
TEMP=$(grep -l coretemp /sys/class/hwmon/hwmon*/name 2>/dev/null | head -1 | xargs -r dirname)
sectors() { awk -v d=$DISK '$3==d{print $6}' /proc/diskstats; }
prev=$(sectors)
while true; do
  sleep 15
  now=$(sectors)
  mbs=$(( (now - prev) * 512 / 15 / 1048576 ))
  prev=$now
  mhz=$(( $(cat /sys/devices/system/cpu/cpu$CORE/cpufreq/scaling_cur_freq) / 1000 ))
  temp=$([ -n "$TEMP" ] && echo $(( $(cat $TEMP/temp1_input) / 1000 )) || echo -)
  read load _ < /proc/loadavg
  mem=$(awk '/MemAvailable/{printf "%.1f", $2/1048576}' /proc/meminfo)
  swap=$(awk '/SwapTotal/{t=$2} /SwapFree/{f=$2} END{printf "%.1f", (t-f)/1048576}' /proc/meminfo)
  echo "$(date +%F_%T) load=$load coreMHz=$mhz tempC=$temp memAvailGB=$mem swapGB=$swap diskMBs=$mbs" >> "$OUT"
done
