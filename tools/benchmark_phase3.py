"""Single-worker Phase 3 public-API timing with deterministic float input.

Run baseline/candidate in separate processes, alternating their order. Each trial
uses a fresh filter, warms its model/FFT/workspace and requests distinct interior
frames. Timings include host/input delivery; output hashing is outside timing.
"""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import sys
import time

import numpy as np
import vapoursynth as vs

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tests/vapoursynth'))
from fixtures import environment


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--plugin', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--trials', type=int, default=3)
    a = p.parse_args()
    if a.trials < 1:
        p.error('--trials must be positive')
    environment(vs)
    c = vs.core
    c.num_threads = 1
    c.std.LoadPlugin(path=str(a.plugin.resolve()))
    samples = np.random.default_rng(142).uniform(0, 1, (16, 360, 640)).astype(np.float32)
    base = c.std.BlankClip(width=640, height=360, length=16, format=vs.GRAYS)

    source_frames = [base.get_frame(n).copy() for n in range(16)]
    for n, frame in enumerate(source_frames):
        np.copyto(np.asarray(frame[0]), samples[n])
    src = c.std.ModifyFrame(base, clips=base, selector=lambda n, f: source_frames[n])
    cases = []
    for algorithm in ('FFT3D', 'DFTTest'):
        for temporal in (1, 3):
            if algorithm == 'FFT3D':
                params = dict(bt=temporal, bw=32, bh=32, ow=16, oh=16,
                              pfactor=1, pframe=0, px=1, py=1, sharpen=.3, dehalo=.2)
            else:
                params = dict(tbsize=temporal, sbsize=16, sosize=12,
                              slocation=[0, 2, 1, 12], ftype=0)
            times, hashes = [], []
            for _ in range(a.trials):
                node = getattr(c.neo_fft, algorithm)(src, opt=0, **params)
                node.get_frame(3)
                start = time.perf_counter()
                frames = [node.get_frame(n) for n in range(4, 12)]
                times.append((time.perf_counter() - start) * 1000 / len(frames))
                hashes.append(hashlib.sha256(b''.join(np.asarray(f[0]).tobytes() for f in frames)).hexdigest())
            assert len(set(hashes)) == 1
            cases.append(dict(algorithm=algorithm, temporal=temporal, params=params,
                              ms_per_frame=times, median_ms=statistics.median(times), output_sha256=hashes[0]))
    info = {k: v.decode() if isinstance(v, bytes) else v for k, v in c.neo_fft.KernelInfo().items()}
    report = dict(plugin_sha256=hashlib.sha256(a.plugin.read_bytes()).hexdigest(),
                  kernel_info=info, source_sha256=hashlib.sha256(samples.tobytes()).hexdigest(),
                  format='GRAYS', width=640, height=360, workers=1, cases=cases)
    a.report.write_text(json.dumps(report, indent=2))
    for case in cases:
        print(case['algorithm'], case['temporal'], f"{case['median_ms']:.3f} ms/frame", flush=True)
    assert len(source_frames) == 16


if __name__ == '__main__':
    main()
