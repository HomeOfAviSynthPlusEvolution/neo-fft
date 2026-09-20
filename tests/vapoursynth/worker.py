import argparse
import hashlib
import json
from pathlib import Path
import platform
import traceback
import numpy as np
import vapoursynth as vs
from fixtures import environment, source

def run(args):
    policy = environment(vs)
    core = vs.core
    core.num_threads = 4
    before = {p.identifier for p in core.plugins()}
    core.std.LoadPlugin(path=str(args.plugin.resolve()))
    loaded = [p for p in core.plugins() if p.identifier not in before]
    if len(loaded) != 1: raise RuntimeError('exactly one explicitly loaded plugin required')
    plugin = loaded[0]
    if Path(plugin.plugin_path).resolve() != args.plugin.resolve(): raise RuntimeError('loaded DLL path differs')
    functions = list(plugin.functions())
    registration = {f.name: f.signature for f in functions}
    if args.new:
        name = args.algorithm
    else:
        candidates = [f.name for f in functions if args.algorithm.lower() in f.name.lower()]
        if len(candidates) != 1: raise RuntimeError(f'ambiguous reference registration: {registration}')
        name = candidates[0]
    call = getattr(plugin, name)
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = dict(plugin=str(args.plugin.resolve()),sha256=hashlib.sha256(args.plugin.read_bytes()).hexdigest(),
        identifier=plugin.identifier,namespace=plugin.namespace,registration=registration,function=name,
        python=platform.python_version(),vs=str(vs.__version__),os=platform.platform(),seed=args.seed,opt=args.opt,
        fft_backend='pocketfft',fft_threads=1,worker_threads=1,autoload=False,cases=[])
    if args.new:
        manifest['kernel_info'] = core.neo_fft.KernelInfo()
        target = manifest['kernel_info']['target']
        manifest['effective_own_kernel'] = 'scalar' if args.opt==1 else target
        if args.opt==0 and target in ('scalar','SCALAR','EMU128','scalar (Highway disabled)'):
            raise RuntimeError('required Highway comparison has no active SIMD target')
    cases = json.loads(args.cases.read_text())
    for case in cases:
        if case['algorithm'] != args.algorithm: continue
        record = dict(case=case)
        try:
            src, inputs, digest = source(vs,case,args.seed)
            record['input_sha256'] = digest
            params = dict(case['params'],fft_backend='pocketfft',opt=args.opt)
            params.update(dict(bt=1,ncpu=1,mt=False) if args.algorithm=='FFT3D' else dict(tbsize=1,threads=1,fft_threads=1))
            if not args.new and params.get('planes') == []:
                try:
                    call(src,**params)
                except vs.Error as e:
                    record['interface_exception'] = str(e)
                else:
                    raise AssertionError('reference empty-array registration changed; review exception')
                if args.algorithm == 'FFT3D':
                    del params['planes']
                    record['reference_mapping'] = 'omitted planes selects all; reference rejects explicit empty list'
                else:
                    record.update(status='interface_exception',reference_mapping='no-plane copy verified against source in candidate; reference registration rejects []')
                    manifest['cases'].append(record)
                    continue
            record['resolved_call'] = params
            dst = call(src,**params)
            if (dst.width,dst.height,dst.num_frames,dst.fps_num,dst.fps_den,dst.format.id) != (src.width,src.height,src.num_frames,src.fps_num,src.fps_den,src.format.id):
                raise AssertionError('output structural metadata changed')
            arrays = {}
            hashes = []
            for n in range(dst.num_frames):
                f = dst.get_frame(n)
                original = src.get_frame(n)
                if dict(f.props) != dict(original.props): raise AssertionError('frame properties changed')
                planes = case['params'].get('planes', list(range(dst.format.num_planes)))
                if not planes and args.algorithm=='FFT3D': planes = list(range(dst.format.num_planes))
                for p in range(dst.format.num_planes):
                    values = np.asarray(f[p]).copy()
                    if not np.isfinite(values).all(): raise AssertionError('non-finite output')
                    if p not in planes and values.tobytes() != inputs[n][p].tobytes(): raise AssertionError('unselected plane changed')
                    arrays[f'{n}-{p}'] = values
                    hashes.append(hashlib.sha256(values.tobytes()).hexdigest())
            if args.new:
                order = sorted(range(dst.num_frames), key=lambda n: (n*5+3)%dst.num_frames)
                order += list(reversed(range(dst.num_frames))) + [0, dst.num_frames-1, 0]
                # A fresh instance has no cached outputs, so distinct frames execute concurrently.
                second = call(src,**params)
                pending = [(n,second.get_frame_async(n)) for n in order]
                for n,future in pending:
                    frame = future.result()
                    for p in range(dst.format.num_planes):
                        if np.asarray(frame[p]).tobytes() != arrays[f'{n}-{p}'].tobytes(): raise AssertionError('request-order mismatch')
            np.savez(args.output/(case['id']+'.npz'),**arrays)
            record.update(status='ok',output_hashes=hashes)
        except Exception:
            record.update(status='error',error=traceback.format_exc())
        manifest['cases'].append(record)
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    failed = [r for r in manifest['cases'] if r['status'] not in ('ok','interface_exception')]
    print(json.dumps(dict(output=str(args.output),cases=len(manifest['cases']),failed=len(failed))))
    for r in failed[:5]: print(r['case']['id'],r['error'])
    return bool(failed)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--plugin',type=Path,required=True); p.add_argument('--cases',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True); p.add_argument('--algorithm',choices=('FFT3D','DFTTest'),required=True)
    p.add_argument('--seed',type=int,required=True); p.add_argument('--opt',type=int,required=True)
    p.add_argument('--new',action='store_true')
    raise SystemExit(run(p.parse_args()))
