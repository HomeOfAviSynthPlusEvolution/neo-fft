"""Phase-3 curve admission, public-path regression and optional black-box comparison.

Frozen limits: integer 1 LSB, float 4e-6 (existing DFTTest budget);
independent rectangular-window center oracle: 4e-6 float samples.
Reference DLL is an opaque input; no reference source inspection is used.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import vapoursynth as vs
from admission import fails
from fixtures import environment, source


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plugin', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--report', type=Path)
    parser.add_argument('--reference-opt', type=int, choices=(0, 1), default=1)
    parser.add_argument('--trace', action='store_true')
    args = parser.parse_args()
    policy = environment(vs)
    core = vs.core
    core.num_threads = 4
    core.std.LoadPlugin(path=str(args.plugin.resolve()))
    call = core.neo_fft.DFTTest
    reference = None
    if args.reference:
        core.std.LoadPlugin(path=str(args.reference.resolve()))
        reference = core.neo_dfttest.DFTTest
    src = core.std.BlankClip(width=32, height=24, format=vs.GRAY8, length=15, color=[64])
    for bad in ([0], [0, 1], [0, 1, .5, 2], [0, -1, 1, 2],
                [0, 1, 0, 2, 1, 3], [0, 1, .5, 2, .500000001, 3, 1, 4],
                [0, 1, 1, float('nan')], [0, 1, 1, 1e100]):
        for key in ('slocation', 'ssx', 'ssy', 'sst'):
            fails(lambda: call(src, planes=[], **{key: bad}))
        fails(lambda: call(src, slocation=[0, 1, 1, 1], ssx=bad))
    fails(lambda: call(src, ssystem=2))
    # Empty planes validates raw input without profile/window construction.
    call(src, planes=[], slocation=[0, 3e38, 1, 3e38]).get_frame(0)
    # Singleton extension: axis-only spatial ordinates do not replace DC.
    for opt in (0, 1):
        single = call(src, tbsize=1, sbsize=1, smode=0, ftype=2, sigma=1,
                      ssx=[0, 99, 1, 99], opt=opt, zmean=False)
        assert np.array_equal(np.asarray(single.get_frame(0)[0]), np.asarray(src.get_frame(0)[0]))

    # Independent NumPy DFT oracle: 3x3x3 rectangular analysis, center synthesis.
    case = dict(format='gray', bits=32, width=24, height=20, frames=5, pattern=6)
    clip, inputs, _ = source(vs, case, 101)
    coeff = np.empty((3, 3, 2))
    for t in range(3):
        for y in range(3):
            for x in range(2):
                r = np.sqrt((min(t, 3-t)**2 + min(y, 3-y)**2 + x*x) / 3)
                coeff[t, y, x] = .2 * (1-r) + .8 * r
    for opt in (0, 1):
        out = call(clip, sbsize=3, smode=0, tbsize=3, swin=7, twin=7, zmean=False,
                   ftype=2, slocation=[0, .2, 1, .8], ssystem=1, opt=opt).get_frame(2)
        actual = np.asarray(out[0])
        for y in range(1, 19):
            for x in range(1, 23):
                block = np.stack([inputs[n][0][y-1:y+2, x-1:x+2] for n in (1, 2, 3)])
                expected = np.fft.irfftn(np.fft.rfftn(block) * coeff, s=(3, 3, 3), axes=(0, 1, 2))[1, 1, 1]
                assert abs(actual[y, x] - expected) <= 4e-6

    records = []
    for bits, family in ((8, 'gray'), (10, '420'), (16, '444'), (32, 'rgb')):
        for T, S in ((1, 8), (3, 5), (5, 8), (15, 3)):
            for pattern in range(7):
                case = dict(format=family, bits=bits, width=32, height=24, frames=max(5, T), pattern=pattern)
                clip, inputs, _ = source(vs, case, 17)
                for ftype in range(5):
                    lo, hi = (2, 12) if ftype < 2 else (.25, .9)
                    variants = (
                        dict(slocation=[1, hi, 0, lo], ssystem=0),
                        dict(slocation=[0, lo, 1, hi], ssystem=1),
                        dict(ssx=[0, lo, 1, hi], ssy=[0, hi, 1, lo], sst=[0, lo, 1, hi]),
                        dict(ssx=[0, lo, 1, hi], ssystem=1),
                    )
                    for variant in variants:
                        params = dict(tbsize=T, sbsize=S, sosize=S//2, ftype=ftype, sigma=lo,
                                      planes=[0], zmean=(ftype % 2 == 0), fft_backend='pocketfft', **variant)
                        if args.trace:
                            print(bits, family, T, S, pattern, ftype, variant, flush=True)
                        outputs = [call(clip, opt=opt, **params) for opt in (1, 0)]
                        # Reference dispatch is recorded separately; candidate always exercises both paths.
                        if reference:
                            outputs.append(reference(clip, opt=args.reference_opt, **params))
                        order = [clip.num_frames-1, 0, clip.num_frames//2]
                        maximum = 0.0
                        for n in order:
                            frames = [out.get_frame(n) for out in outputs]
                            arrays = [np.asarray(f[0]).astype(float) for f in frames]
                            limit = 4e-6 if bits == 32 else 1
                            for a in arrays[1:]:
                                delta = float(np.max(np.abs(a-arrays[0])))
                                maximum = max(maximum, delta)
                                assert delta <= limit, (bits, family, T, S, ftype, variant, n, delta)
                            for frame in frames:
                                assert dict(frame.props) == dict(clip.get_frame(n).props)
                                for p in range(1, clip.format.num_planes):
                                    assert np.asarray(frame[p]).tobytes() == inputs[n][p].tobytes()
                        # Fresh candidate instances, concurrent first calls and repeatability.
                        fresh = call(clip, opt=0, **params)
                        pending = [(n, fresh.get_frame_async(n)) for n in order]
                        for n, future in pending:
                            assert np.array_equal(np.asarray(future.result()[0]), np.asarray(outputs[1].get_frame(n)[0]))
                        records.append(dict(pattern=pattern, bits=bits, family=family, T=T, S=S, ftype=ftype,
                                            curves=variant, max_error=maximum))
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    info = {k: v.decode() if isinstance(v, bytes) else v for k, v in core.neo_fft.KernelInfo().items()}
    report = dict(candidate_sha256=digest(args.plugin), reference_sha256=digest(args.reference) if args.reference else None,
                  python=platform.python_version(), vs=str(vs.__version__), kernel_info=info,
                  reference_opt=args.reference_opt, reference_dispatch='not instrumented', cases=records)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2))
    print(f'Phase-3 curves: admission, independent DFT, {len(records)} public cases passed; {info}')


if __name__ == '__main__':
    main()
