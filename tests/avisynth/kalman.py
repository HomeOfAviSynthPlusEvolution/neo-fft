"""AVS warmup parameter and cold-seek dependency bounds with remapped oracles."""


def check_kalman(run,opt):
    for bits in (8,16,32):
        setup=('c=ColorBars(width=128,height=96).ConvertToYV12().ExtractY().Trim(0,7)\n'
               f'c=(c+c.Invert()).Loop(1000).ConvertBits({bits})\n')
        for warmup in (0,4,8,16):
            # Full-prefix reference local frame zero maps to just before warmup;
            # cold requested local W+1 maps to original frame 10000.
            options=f'bt=0,bw=8,bh=8,sigma=12,opt={opt}'
            body=f'''guard=c.ScriptClip("""Assert(current_frame>={10000-warmup} && current_frame<=10002,"outside warmup budget")
last""")
global expected=neo_fft_FFT3D(c.Trim({9999-warmup},0),{options},kalman_warmup=2147483647).Trim({warmup+1},0)
o=neo_fft_FFT3D(guard,{options},kalman_warmup={warmup}).Trim(10000,0)
return o.ScriptClip("""Assert(LumaDifference(last,expected)==0,"warmup differs from bounded prefix")
last""")
'''
            run(f'kalman-warmup-{bits}-{warmup}',body,prefix=setup,frame=0)
        # Default is eight, and the newly appended parameter accepts Undefined().
        body=f'''global expected=neo_fft_FFT3D(c,bt=0,bw=8,bh=8,opt={opt},kalman_warmup=8)
o=neo_fft_FFT3D(c,bt=0,bw=8,bh=8,opt={opt},kalman_warmup=Undefined())
return o.ScriptClip("""Assert(LumaDifference(last,expected)==0,"default warmup differs")
last""")
'''
        run(f'kalman-default-{bits}',body,prefix=setup,frame=10000)
    for value in (-1,2147483648):
        run(f'kalman-invalid-{value}',f'return neo_fft_FFT3D(c,bt=1,planes=[],kalman_warmup={value})',error='kalman_warmup')
