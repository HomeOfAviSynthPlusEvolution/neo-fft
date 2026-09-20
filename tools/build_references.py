"""Build the fixed Windows reference DLLs using GCC; never changes reference kernels.

SDKs and toolchain are explicit arguments. An isolated header overlay restores
the exact historical DS2 mdspan API still used by the pinned FFT3D source.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

FFT3D='bc9dd39410d06470d55ffd9147a98b9252f451f3'
DFTTEST='e5d59aaff79aa45a8dc7ccb5501aff39dd489475'
DS2='2a24d6b4fe808692bfa10c1f9734a3c50c15774e'
MDSPAN='e0a68c580a967ffc017708e35bab5514b71daf42'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--dualsynth',type=Path,required=True)
    p.add_argument('--pocketfft',type=Path,required=True)
    p.add_argument('--gcc-bin',type=Path,required=True)
    p.add_argument('--vs-include',type=Path,required=True,help='contains vapoursynth/VapourSynth4.h')
    p.add_argument('--avs-include',type=Path,required=True,help='contains avisynth.h and avisynth_c.h')
    args=p.parse_args(); root=args.root.resolve(); root.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ); env['PATH']=str(args.gcc_bin.resolve())+os.pathsep+env['PATH']
    recipe=[]
    def run(cmd):
        cmd=[str(x) for x in cmd];recipe.append(cmd)
        subprocess.run(cmd,env=env,check=True)
    def git(*argv): return subprocess.check_output(['git',*map(str,argv)],env=env)
    if git('-C',args.dualsynth,'rev-parse','HEAD').decode().strip()!=DS2: raise RuntimeError('DS2 revision differs from spec')
    if git('-C',args.dualsynth,'diff','--name-only','HEAD','--','include','src').strip(): raise RuntimeError('DS2 implementation is dirty')
    if git('-C',args.pocketfft,'rev-parse','HEAD').decode().strip()!='5f27d5a8f51c5c25030cb22abf434decc9faf0ff':
        raise RuntimeError('PocketFFT revision differs')
    overlay=root/'compat/dualsynth/mdspan.hpp';overlay.parent.mkdir(parents=True,exist_ok=True)
    overlay.write_bytes(git('-C',args.dualsynth,'show',MDSPAN+':include/dualsynth/mdspan.hpp'))
    for name,repo,commit,target in [('fft3d','neo_FFT3D',FFT3D,'neo-fft3d'),('dfttest','neo_DFTTest',DFTTEST,'neo-dfttest')]:
        src=root/(name+'-src');build=root/(name+'-build')
        if not src.exists(): run(['git','clone','https://github.com/HomeOfAviSynthPlusEvolution/'+repo+'.git',src])
        run(['git','-C',src,'checkout','--detach',commit])
        flags=f'-I"{args.vs_include.resolve().as_posix()}" -I"{args.avs_include.resolve().as_posix()}"'
        if name=='fft3d': flags+=f' -I"{(root/"compat").as_posix()}" -DDS_USE_STD_MDSPAN=1 -include windows.h -DNOMINMAX'
        run(['cmake','-S',src,'-B',build,'-G','Ninja','-DCMAKE_BUILD_TYPE=Release',
             '-DCMAKE_CXX_COMPILER='+str((args.gcc_bin/'g++.exe').resolve()),
             '-DFETCHCONTENT_SOURCE_DIR_DUALSYNTH2='+str(args.dualsynth.resolve()),
             '-DFETCHCONTENT_SOURCE_DIR_POCKETFFT='+str(args.pocketfft.resolve()),
             '-DCMAKE_SHARED_LINKER_FLAGS=-static','-DCMAKE_CXX_FLAGS='+flags])
        run(['cmake','--build',build,'--parallel','4','--target',target])
    manifest=dict(recipe=recipe,ds2=DS2,fft3d=FFT3D,dfttest=DFTTEST,mdspan_source=MDSPAN,
                  overlay_sha256=hashlib.sha256(overlay.read_bytes()).hexdigest(),binaries={})
    for file in root.glob('*-build/*neo*.dll'):
        manifest['binaries'][str(file)]=hashlib.sha256(file.read_bytes()).hexdigest()
    (root/'build-manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
