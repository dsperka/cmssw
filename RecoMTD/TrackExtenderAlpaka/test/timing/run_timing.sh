#!/bin/bash
# Timing of the full reconstruction: legacy MTD timing, Alpaka MTD timing, Alpaka MTD timing + device 3D
# vertices and TOF PID; plus a read-only job that gives the most the input can deliver.
# Run inside cmsenv from an empty work directory:
#
#   run_timing.sh [maxEvents] [threads] [fileList] [filePrefix] [cores]
#
#   maxEvents   default 300
#   threads     threads = streams per job, default 1
#   fileList    GEN-SIM-RECO LFNs, default ../validation/relvalTTbar_PU200_D128_20_1_0_pre3.txt
#   filePrefix  prepended to each LFN, default root://cms-xrd-global.cern.ch/ (use file:/path/to/store/..
#               for local copies: the page cache can then be dropped before each job)
#   cores       CPU list for taskset, default 2 (one physical core; leave its hyperthread sibling idle)
#
# The jobs run one after the other, never at the same time. Before each job the local input files are
# dropped from the page cache and the script waits until the machine is idle (1-minute load < 1, then
# one more minute). monitor.sh logs the machine state every 15 s to monitor.log.
# Output per job: <job>.log (ThroughputService), <job>.json (FastTimerService, per module), fjr_<job>.xml.
# Summary: python3 timing_summary.py .
set -e
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
NEV=${1:-300}
THREADS=${2:-1}
LIST=${3:-$HERE/../validation/relvalTTbar_PU200_D128_20_1_0_pre3.txt}
PREFIX=${4:-root://cms-xrd-global.cern.ch/}
CORES=${5:-2}
FILES=$(sed "s|^|${PREFIX}|" "$LIST" | paste -sd, -)
LOCAL=""
[[ $PREFIX == file:* ]] && LOCAL=$(sed "s|^|${PREFIX#file:}|" "$LIST")

CORE=${CORES%%[,-]*} "$HERE/monitor.sh" monitor.log &
MON=$!
trap 'kill $MON' EXIT

settle() {
  [ -n "$LOCAL" ] && python3 "$HERE/drop_page_cache.py" $LOCAL
  until awk '{exit !($1 < 1.0)}' /proc/loadavg; do sleep 30; done
  sleep 60
  echo "=== $1 start $(date +%F_%T) load=$(cut -d' ' -f1-3 /proc/loadavg)" | tee -a monitor.log
}

# the input branches the reconstruction reads, from a short job's framework job report
if [ ! -s branches.txt ]; then
  settle probe
  taskset -c $CORES cmsRun -j fjr_probe.xml "$HERE/mtdTiming_cfg.py" timing=legacy threads=1 maxEvents=3 \
    inputFiles="$FILES" output=probe.root > probe.log 2>&1
  python3 - <<'EOF'
import xml.etree.ElementTree as ET
names = sorted(b.get('Name') for rb in ET.parse('fjr_probe.xml').getroot().iter('ReadBranches') for b in rb)
open('branches.txt', 'w').write('\n'.join(names) + '\n')
EOF
fi

settle readonly
taskset -c $CORES cmsRun -j fjr_readonly.xml "$HERE/readOnly_cfg.py" inputFiles="$FILES" maxEvents=$NEV \
  branches=branches.txt > readonly.log 2>&1
echo "readonly done $(date +%F_%T)"

for job in "legacy:timing=legacy" "alpaka:timing=alpaka" "device:timing=alpaka vertices=alpaka"; do
  name=${job%%:*}
  settle $name
  taskset -c $CORES cmsRun -j fjr_$name.xml "$HERE/mtdTiming_cfg.py" ${job#*:} threads=$THREADS maxEvents=$NEV \
    inputFiles="$FILES" output=$name.root > $name.log 2>&1
  echo "$name done $(date +%F_%T)"
done
python3 "$HERE/timing_summary.py" . | tee summary.txt
