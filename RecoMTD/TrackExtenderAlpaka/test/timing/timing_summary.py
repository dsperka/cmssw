#!/usr/bin/env python3
"""Tables from the output of run_timing.sh: throughput, time of the modules the MTD chain replaces,
time of all other modules (a check of the machine stability), and the machine state during each job.

usage: python3 timing_summary.py workdir
"""
import json
import os
import re
import statistics
import sys

JOBS = ['legacy', 'alpaka', 'device']
NAMES = {'readonly': 'Read only', 'legacy': 'Legacy', 'alpaka': 'Alpaka timing',
         'device': 'Alpaka timing + device vertices, TOF PID'}
# the modules of the MTD timing, the 3D vertices and tofPID4DnoPID, legacy and Alpaka
STEPS = [('MTD timing + quality BDT', ['trackExtenderWithMTD', 'mtdTrackQualityMVA', 'mtdEventInput', 'mtdTimingProducer']),
         ('3D vertices (both fits)', ['unsortedOfflinePrimaryVertices', 'pvVertexSoA',
                                      'unsortedOfflinePrimaryVerticesFromSoA', 'unsortedOfflinePrimaryVerticesFromSoAWithBS']),
         ('TOF PID (tofPID4DnoPID)', ['tofPID4DnoPID', 'tofPidProducer'])]
NOT_MODULES = {'other', 'eventsetup', 'source'}


def throughput(log):
    for line in open(log):
        m = re.match(r'Average throughput: ([0-9.]+) ± ([0-9.]+) ev/s', line)
        if m:
            return float(m.group(1)), float(m.group(2))
    return None


def machine_state(path):
    state, job = {}, None
    for line in open(path):
        m = re.match(r'=== (\w+) start', line)
        if m:
            job = m.group(1)
            state[job] = []
        elif job and 'coreMHz=' in line:
            v = dict(f.split('=', 1) for f in line.split()[1:] if '=' in f)
            state[job].append((float(v['load']), int(v['coreMHz']), v['tempC']))
    return state


def main():
    os.chdir(sys.argv[1] if len(sys.argv) > 1 else '.')
    print('Throughput (ThroughputService, first events skipped)')
    for job in ['readonly'] + JOBS:
        if os.path.exists(job + '.log') and throughput(job + '.log'):
            t, e = throughput(job + '.log')
            print('  %-42s %.5f ± %.5f ev/s   %7.2f s/event' % (NAMES[job], t, e, 1. / t))

    timing = {j: json.load(open(j + '.json')) for j in JOBS if os.path.exists(j + '.json')}
    per_event = {j: {m['label']: m['time_thread'] / d['total']['events'] for m in d['modules']}
                 for j, d in timing.items()}
    print('\nTime per event of the replaced modules (FastTimerService, thread time, ms)')
    print('  %-28s' % '' + ''.join('%14s' % j for j in timing))
    chain = set()
    total = {j: 0. for j in timing}
    for step, labels in STEPS:
        chain.update(labels)
        row = {j: sum(per_event[j].get(label, 0.) for label in labels) for j in timing}
        for j in timing:
            total[j] += row[j]
        print('  %-28s' % step + ''.join('%14.1f' % row[j] for j in timing))
    print('  %-28s' % 'total' + ''.join('%14.1f' % total[j] for j in timing))
    print('  %-28s' % 'all other modules' + ''.join(
        '%14.1f' % sum(t for label, t in per_event[j].items() if label not in chain | NOT_MODULES) for j in timing))

    if os.path.exists('monitor.log'):
        print('\nMachine state during each job (monitor.log, every 15 s)')
        for job, samples in machine_state('monitor.log').items():
            if not samples:
                continue
            mhz = [s[1] for s in samples]
            load = [s[0] for s in samples]
            temps = [int(s[2]) for s in samples if s[2] != '-']
            print('  %-10s core MHz mean %4.0f median %4.0f, below 3 GHz %4.1f%% | load mean %.2f max %.2f | temp max %s C'
                  % (job, statistics.mean(mhz), statistics.median(mhz), 100. * sum(f < 3000 for f in mhz) / len(mhz),
                     statistics.mean(load), max(load), max(temps) if temps else '-'))


if __name__ == '__main__':
    main()
