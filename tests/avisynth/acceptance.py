"""Focused host tests for the AviSynth C++ adapter, using the pinned DS2 runner."""
import argparse
import ctypes
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    for name in ("runner", "plugin", "runtime", "work"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--host-api", choices=("c", "cpp"), default="c")
    parser.add_argument("--opt", type=int, choices=(0, 1), default=1)
    parser.add_argument("--planes-only", action="store_true")
    parser.add_argument("--kalman-only", action="store_true")
    args = parser.parse_args()
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    plugin = Path(args.plugin).resolve().as_posix()
    if os.name == "nt":
        module = ctypes.WinDLL(plugin)
        assert hasattr(module, "AvisynthPluginInit3") or hasattr(module, "_AvisynthPluginInit3@8")
        assert not hasattr(module, "avisynth_c_plugin_init2")
    header = f'LoadPlugin("{plugin}")\n'
    fixture = 'c=BlankClip(width=128,height=96,length=8,fps=24,pixel_type="YV12",color_yuv=$608080)\n'
    count = 0

    def run(name, body, *, prefix=fixture, frame=4, error=None, mode="--video", extra=()):
        nonlocal count
        script = work / (name + ".avs")
        script.write_text(header + prefix + body, encoding="utf-8")
        result = subprocess.run(
            [args.runner, mode, str(script), "--backend", args.host_api, "--runtime", args.runtime,
             "--frame", str(frame), *extra], capture_output=True, text=True, timeout=40)
        text = result.stdout + result.stderr
        assert (result.returncode != 0 and error in text) if error else result.returncode == 0, (name, text)
        count += 1
        return result.stdout

    opt = args.opt
    if args.kalman_only:
        from kalman import check_kalman
        check_kalman(run,opt)
        print(f"AviSynth Kalman warmup: {count} cases passed (opt={opt})")
        return
    if args.planes_only:
        from planes import check_planes
        check_planes(run, opt)
        print(f"AviSynth plane selection: {count} cases passed (opt={opt})")
        return
    run("registration", '''
Assert(FunctionExists("neo_fft_FFT3D"))
Assert(FunctionExists("neo_fft_DFTTest"))
Assert(FunctionExists("neo_fft_KernelInfo"))
k=neo_fft_KernelInfo()
Assert(IsString(k[0]) && IsString(k[1]) && IsString(k[2]) && IsInt(k[3]))
Assert(k[3]>=1)
return c
''')
    for bt in range(-1, 6):
        run(f"fft-bt{bt}", f'return neo_fft_FFT3D(c,bt={bt},opt={opt}).Prefetch(4)')
    variants = {
        "fft-arrays-roi": f'neo_fft_FFT3D(c,planes=[0,2],l=2,r=2,t=2,b=2,interlaced=true,opt={opt})',
        "fft-scalar-plane": f'neo_fft_FFT3D(c,planes=0,cache_frames=2,cache_mb=1,opt={opt})',
        "fft-preview": f'neo_fft_FFT3D(c,pfactor=1,pshow=true,px=1,py=1,opt={opt})',
        "fft-sampled": f'neo_fft_FFT3D(c,pfactor=1,pframe=7,px=1,py=1,opt={opt})',
        "dft-default": f'neo_fft_DFTTest(c,opt={opt})',
        "dft-ola": f'neo_fft_DFTTest(c,tmode=1,tbsize=4,tosize=2,opt={opt})',
        "dft-curves": f'neo_fft_DFTTest(c,slocation=[0,4,1,12],ssx=[0,1,1,2],ssystem=1,opt={opt})',
        "dft-sampled": f'neo_fft_DFTTest(c,nlocation=[0,0,0,0],alpha=5,opt={opt})',
        "dft-dither": f'neo_fft_DFTTest(c,planes=0,dither=2,dither_seed=7,threads=16,opt={opt})',
        "dft-empty": f'neo_fft_DFTTest(c,planes=[],opt={opt})',
    }
    for name, expression in variants.items():
        run(name, "return " + expression + ".Prefetch(4)")
    # String parsing is an AVS adapter feature; both forms must feed identical
    # numeric arrays to the shared core, on spatially and temporally varying data.
    string_pairs = {
        "shared": ('slocation=[0,0.25,1,0.75]', 'slocation="0 +.25 1 7.5e-1"'),
        "comma": ('slocation=[0,0.25,1,0.75]', 'slocation="0,+.25,1,7.5e-1"'),
        "separator-axes": ('ssx=[0,0.5,1,2],ssy=[0,1,1,3],sst=[0,2,1,4]',
                           'ssx="0:.5 1:2",ssy=",0,1:1,3:",sst="0::2,, 1:4"'),
        "separator-sample": ('nlocation=[0,0,0,0,2,1,2,4]', 'nlocation=":+0,0:0,0 ,: 2,1:2,4,"'),
        "axes": ('ssx=[0,0.5,1,2],ssy=[0,1,1,3],sst=[0,2,1,4]',
                 'ssx="0 .5 1 2",ssy="0 1 1 3",sst="0 2 1 4"'),
        "mixed": ('ssx=[0,0.5,1,2],ssy=[0,1,1,3]', 'ssx="0 .5 1 2",ssy=[0,1,1,3]'),
        "sample": ('nlocation=[0,0,0,0,2,1,2,4]', 'nlocation="+0 0 0 0 2 1 2 4"'),
        "empty": ('nlocation=[],slocation=[],ssx=[],ssy=[],sst=[]',
                  'nlocation=" ,:",slocation=" : , ",ssx="",ssy=",,",sst="::"'),
        "whitespace": ('slocation=[0,0.25,1,0.75]',
                       'slocation=" 0"+Chr(9)+".25"+Chr(13)+Chr(10)+"1 .75 "'),
    }
    for pixel, convert in (("u8", ""), ("u16", ".ConvertBits(16)"), ("float", ".ConvertBits(32)")):
        setup = ('c=ColorBars(width=128,height=96).ConvertToYV12().Trim(0,3)\n'
                 'c=(c+c.Invert())' + convert + '\n')
        for name, (array_args, text_args) in string_pairs.items():
            body = f'''global expected=neo_fft_DFTTest(c,{array_args},opt={opt})
o=neo_fft_DFTTest(c,{text_args},opt={opt})
return o.ScriptClip("""Assert(LumaDifference(last,expected)==0)
Assert(ChromaUDifference(last,expected)==0 && ChromaVDifference(last,expected)==0)
last""").Prefetch(4)
'''
            for frame in (0, 4, 7):
                run(f"string-{pixel}-{name}-{frame}", body, prefix=setup, frame=frame)
    for pixel in ("Y8", "YUV420P10", "YUV420P16", "Y32", "YUV444PS", "RGBP", "RGBPS"):
        setup = f'c=BlankClip(width=128,height=96,length=8,pixel_type="{pixel}")\n'
        for function in ("FFT3D", "DFTTest"):
            run(f"{function}-{pixel}", f'return neo_fft_{function}(c,opt={opt})', prefix=setup)

    # A temporal volume mean catches slot/center mistakes on nonconstant input.
    temporal = '''
c=BlankClip(width=16,height=16,length=1,pixel_type="Y8",color_yuv=$018080) \
 + BlankClip(width=16,height=16,length=1,pixel_type="Y8",color_yuv=$028080) \
 + BlankClip(width=16,height=16,length=1,pixel_type="Y8",color_yuv=$068080)
'''
    mean = f'neo_fft_DFTTest(c,sbsize=1,smode=0,tbsize=3,swin=7,twin=7,ftype=2,sigma=0,opt={opt})'
    run("temporal-mean", "return " + mean, prefix=temporal, frame=1,
        extra=("--expect-y8-sum", str(16 * 16 * 3)))
    # Repeated Kalman requests through one instance retain the same luma after a backward seek.
    # Use built-in per-frame constants; no external source or grain plugin is required.
    sequence = 'c=' + ' + '.join(
        f'BlankClip(width=32,height=32,length=1,pixel_type="Y8",color_yuv=${v:02X}8080)'
        for v in (10, 30, 90, 20, 70, 40, 60, 80)) + '\n'
    kalman = f'o=neo_fft_FFT3D(c,bw=8,bh=8,bt=0,opt={opt})\n'
    run("kalman-late", kalman + "return o", prefix=sequence, frame=6)
    run("kalman-reordered", kalman + '''
return o.ScriptClip("""a=AverageLuma(last,offset=2)
b=AverageLuma(last,offset=-2)
c=AverageLuma(last,offset=2)
Assert(a==c,"Kalman repeated request changed")
last""")
''', prefix=sequence, frame=4)

    for function in ("FFT3D", "DFTTest"):
        expression = f'neo_fft_{function}(c,opt={opt})'
        audio = fixture + 'c=AudioDub(c,Tone(length=1.0,samplerate=48000,channels=2)).AssumeTFF()\n'
        body = f'''o={expression}
Assert(AudioRate(o)==AudioRate(c) && AudioChannels(o)==AudioChannels(c))
Assert(GetParity(o,4)==GetParity(c,4))
'''
        source = run(function + "-audio-source", body + "return c", prefix=audio, mode="--audio")
        output = run(function + "-audio-output", body + "return o.Prefetch(4)", prefix=audio, mode="--audio")
        assert source == output, (function, source, output)
        run(function + "-props", f'''c=c.propSet("NeoFFTTest",123)
o={expression}
return o.ScriptClip("""Assert(propGetInt(last,"NeoFFTTest")==123)
last""").Prefetch(4)
''')

    invalid = [
        ('neo_fft_FFT3D(c,planes=[3])', 'planes'),
        ('neo_fft_DFTTest(c,planes=[0.5])', 'integer'),
        ('neo_fft_DFTTest(c,slocation=[0,1,0,2])', 'curve'),
        ('neo_fft_DFTTest(c,tbsize=4)', 'odd'),
        ('neo_fft_FFT3D(c,cache_mb=-2)', 'cache'),
        ('neo_fft_DFTTest(c,nlocation="0 0 0 1.5")', 'nlocation'),
        ('neo_fft_DFTTest(c,nlocation="0 0 0 2147483648")', 'nlocation'),
        ('neo_fft_DFTTest(c,nlocation="0 0 0")', 'quadruples'),
        ('neo_fft_DFTTest(c,nlocation="0 0 0 -1")', 'negative sample'),
        ('neo_fft_DFTTest(c,slocation="0 1 1")', 'pairs'),
        ('neo_fft_DFTTest(c,slocation="0 1 0 2")', 'curve'),
        ('neo_fft_DFTTest(c,slocation="0 1 1 1e39")', 'float'),
        ('neo_fft_DFTTest(c,slocation="0 1 1 NaN")', 'slocation'),
        ('neo_fft_DFTTest(c,slocation="0 1 1 inf")', 'slocation'),
        ('neo_fft_DFTTest(c,slocation="0 1 1 0x1p2")', 'slocation'),
        ('neo_fft_DFTTest(c,ssx="0 1 1 2junk")', 'ssx'),
        ('neo_fft_DFTTest(c,ssy="0;1,1,2")', 'ssy'),
        ('neo_fft_DFTTest(c,sst="0:1 1:2junk")', 'sst'),
        ('neo_fft_DFTTest(c,slocation="0 1 1 2",ssx="bad")', 'ssx'),
        ('neo_fft_DFTTest(c,planes="0")', 'integer'),
    ]
    for index, (expression, error) in enumerate(invalid):
        run(f"invalid-{index}", "return " + expression, error=error)
    setup = ('c=ColorBars(width=128,height=96).ConvertToYV12().Trim(0,3)\n'
             'c=c+c.Invert()\n')
    for function, ignored in (
        ("FFT3D", 'mt=true,ncpu=2147483647,measure=true,fft_backend="gpu"'),
        ("FFT3D", 'mt=false,ncpu=-2147483647,measure=false,fft_backend="unknown"'),
        ("DFTTest", 'threads=2147483647,fft_threads=-2147483647,fft_backend="fftw"'),
        ("DFTTest", 'threads=-2147483647,fft_threads=2147483647,fft_backend="unknown"')):
        body = f'''global expected=neo_fft_{function}(c,opt={opt})
o=neo_fft_{function}(c,{ignored},opt={opt})
return o.ScriptClip("""Assert(LumaDifference(last,expected)==0)
Assert(ChromaUDifference(last,expected)==0 && ChromaVDifference(last,expected)==0)
last""").Prefetch(4)
'''
        run(function + "-ignored-" + str(count), body, prefix=setup)
    print(f"AviSynth opt={opt}: {count} host cases passed")


if __name__ == "__main__":
    main()
