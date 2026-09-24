"""Plane selection and full-resolution Alpha against independent gray oracles."""


def check_planes(run, opt):
    def compare(name, setup, actual, expected, channels="YUVA", frame=4, reordered=False):
        checks = "\n".join(
            f'Assert(LumaDifference(last.Extract{p}(),expected.Extract{p}())==0,"plane {p}")'
            for p in channels)
        if reordered:
            checks += '''
a=AverageLuma(last.ExtractA(),offset=2)
b=AverageLuma(last.ExtractA(),offset=-2)
c=AverageLuma(last.ExtractA(),offset=2)
Assert(a==c,"Alpha Kalman replay changed")'''
        run(name, f'''global expected={expected}
o={actual}
return o.ScriptClip("""{checks}
last""").Prefetch(4)
''', prefix=setup, frame=frame)

    for family, convert, channels in (("420", ".ConvertToYV12()", "YUVA"),
                                      ("444", ".ConvertToYV24()", "YUVA"),
                                      ("rgb", ".ConvertToPlanarRGB()", "RGBA")):
        for bits in (8, 16, 32):
            setup = ('c=ColorBars(width=128,height=96)' + convert + '.Trim(0,3)\n'
                     f'c=(c+c.Invert()).ConvertBits({bits})\n'
                     f'c=c.AddAlphaPlane(c.Extract{channels[0]}().Invert())\n')
            for function, options in (("FFT3D", "bw=16,bh=16,bt=3,sigma=32,degrid=0"),
                                      ("DFTTest", "sbsize=8,sosize=4,ftype=2,sigma=0.5,zmean=false")):
                call = f"neo_fft_{function}"
                common = f"{options},opt={opt}"
                tag = f"{function}-{family}-{bits}"
                compare(tag + "-default", setup, f"{call}(c,{common})",
                        f"{call}(c,{common},planes=[0,1,2])", channels)
                compare(tag + "-alpha", setup, f"{call}(c,{common},y=2,u=2,v=2,a=3)",
                        f"c.RemoveAlphaPlane().AddAlphaPlane({call}(c.ExtractA(),{common}))", channels)
                compare(tag + "-copy", setup, f"{call}(c,{common},y=2,u=1,v=2,a=1)", "c", channels)

    setup = ('c=ColorBars(width=128,height=96).ConvertToYV12().Trim(0,3)\n'
             'c=c+c.Invert()\nc=c.AddAlphaPlane(c.ExtractY().Invert())\n')
    for function in ("FFT3D", "DFTTest"):
        call = f"neo_fft_{function}"
        for name, actual, expected in (
            ("missing-modes", "y=2,a=3", "planes=[1,2,3]"),
            ("missing-alpha", "y=3", "planes=[0,1,2]"),
            ("precedence", "planes=[3],y=0,u=9,v=-1,a=2", "planes=[3]"),
            ("empty", "planes=[],y=2,u=2,v=2,a=3", "planes=[]"),
        ):
            compare(function + "-" + name, setup, f"{call}(c,{actual},opt={opt})",
                    f"{call}(c,{expected},opt={opt})")
        compare(function + "-allcopy", setup, f"{call}(c,y=2,u=2,v=2,a=2,opt={opt})", "c")
        ignored = ('mt=true,ncpu=-7,measure=true,fft_backend="unused"' if function == "FFT3D"
                   else 'threads=-7,fft_threads=64,fft_backend="unused",ssx="0 2 1 8"')
        shared = "" if function == "FFT3D" else ",ssx=[0,2,1,8]"
        compare(function + "-legacy-with-compat", setup,
                f"{call}(c,y=2,u=2,v=2,a=3,{ignored},opt={opt})",
                f"{call}(c,planes=[3]{shared},opt={opt})")
        # Gray and three-plane clips ignore valid selectors for missing planes.
        for source, channels in (("c.RemoveAlphaPlane()", "YUV"), ("c.ExtractY()", "Y")):
            compare(function + "-absent-" + channels, setup,
                    f"{call}({source},y=2,u=2,v=2,a=3,opt={opt})", source, channels)
        run(function + "-bad-mode", f"return {call}(c,a=4,opt={opt})", prefix=setup,
            error="a: mode must be 1, 2 or 3")
        run(function + "-bad-plane", f"return {call}(c,planes=[4],opt={opt})", prefix=setup,
            error="planes index outside actual format")

    for tag, options in (("kalman", "bt=0,bw=16,bh=16"),
                         ("pattern", "bt=3,bw=16,bh=16,pfactor=1,pframe=2,px=1,py=1"),
                         ("roi", "bt=3,bw=16,bh=16,interlaced=true,l=2,r=2,t=2,b=2")):
        for frame in (0, 4, 7):
            compare("fft-alpha-" + tag + str(frame), setup,
                    f"neo_fft_FFT3D(c,{options},planes=[3],opt={opt})",
                    f"c.RemoveAlphaPlane().AddAlphaPlane(neo_fft_FFT3D(c.ExtractA(),{options},opt={opt}))",
                    frame=frame, reordered=tag == "kalman" and frame == 4)
    compare("fft-alpha-cache-off", setup, f"neo_fft_FFT3D(c,planes=[3],bt=5,cache_mb=0,opt={opt})",
            f"neo_fft_FFT3D(c,planes=[3],bt=5,opt={opt})")
    # The sample rectangle fits Alpha but exceeds the subsampled chroma bounds.
    for temporal in ("tbsize=3", "tbsize=4,tmode=1,tosize=2"):
        options = f"sbsize=8,sosize=4,{temporal},opt={opt}"
        compare("dft-alpha-sample-" + temporal[:8], setup,
                f"neo_fft_DFTTest(c,{options},planes=[3],nlocation=[0,3,60,80])",
                f"c.RemoveAlphaPlane().AddAlphaPlane(neo_fft_DFTTest(c.ExtractA(),{options},nlocation=[0,0,60,80]))")
    compare("dft-alpha-dither", setup,
            f"neo_fft_DFTTest(c,y=2,u=2,v=2,a=3,dither=2,dither_seed=7,opt={opt})",
            f"neo_fft_DFTTest(c,planes=[3],dither=2,dither_seed=7,opt={opt})")
