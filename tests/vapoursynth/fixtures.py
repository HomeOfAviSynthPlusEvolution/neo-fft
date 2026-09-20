"""P1 fixture v1: frame, plane, row, column traversal; uint32 LCG."""
import hashlib
import threading
import numpy as np

def environment(vs):
    class Policy(vs.EnvironmentPolicy):
        def on_policy_registered(self, api):
            self.api = api
            self.env = api.create_environment(int(vs.CoreCreationFlags.DISABLE_AUTO_LOADING))
            self.local = threading.local()
        def get_current_environment(self):
            return getattr(self.local, 'env', self.env)
        def set_environment(self, env):
            previous = self.get_current_environment()
            self.local.env = env
            return previous
        def on_policy_cleared(self):
            self.api.destroy_environment(self.env)
    policy = Policy()
    vs.register_policy(policy)
    return policy

def source(vs, case, seed):
    family, bits = case['format'], case['bits']
    color = vs.GRAY if family == 'gray' else vs.RGB if family == 'rgb' else vs.YUV
    sw, sh = (1, 1) if family == '420' else (1, 0) if family == '422' else (0, 0)
    fmt = vs.core.query_video_format(color, vs.FLOAT if bits == 32 else vs.INTEGER, bits, sw, sh)
    w, h = case.get('width', 256 if sw else 128), case.get('height', 192 if sw else 96)
    count = case.get('frames', 1)
    base = vs.core.std.BlankClip(width=w, height=h, format=fmt.id, length=count, fpsnum=24000, fpsden=1001)
    state = seed
    frames = []
    for n in range(count):
        planes = []
        for p in range(fmt.num_planes):
            ph, pw = h >> (sh if p else 0), w >> (sw if p else 0)
            noise = np.empty(ph * pw, dtype=np.uint32)
            for i in range(noise.size):
                state = (1664525 * state + 1013904223) & 0xffffffff
                noise[i] = state
            yy, xx = np.indices((ph, pw))
            pattern = case.get('pattern', n if count == 7 else 6)
            if pattern == 0: v = np.zeros((ph, pw))
            elif pattern == 1: v = np.full((ph, pw), 0.5 if p else 0.25)
            elif pattern == 2:
                v = np.zeros((ph, pw)); v[0, 0] = 1; v[-1, -1] = 1; v[ph//2, pw//2] = 1
            elif pattern == 3: v = (xx / max(pw-1, 1) + yy / max(ph-1, 1)) / 2
            elif pattern == 4: v = ((xx + yy) % 2).astype(float)
            elif pattern == 5: v = .5 + .2*np.sin(2*np.pi*xx/9) + .2*np.cos(2*np.pi*yy/7)
            elif pattern == 7: v = np.take(np.array([-.25, 0, .5, 1, 1.25]), (xx+yy) % 5)
            elif pattern == 8: v = xx / max(pw-1, 1)
            elif pattern == 9: v = yy / max(ph-1, 1)
            elif pattern == 10: v = .5 + .2*np.sin(2*np.pi*xx/9)
            elif pattern == 11: v = .5 + .2*np.cos(2*np.pi*yy/7)
            else: v = (noise.reshape(ph, pw) >> 8).astype(np.float64) / (1 << 24)
            if bits == 32:
                if family in ('444', '422', '420') and p: v = v - .5
                if pattern == 6: v = v*1.5 - .25
                arr = v.astype(np.float32)
            else:
                arr = np.floor(v*((1 << bits)-1)+.5).astype(np.uint8 if bits == 8 else np.uint16)
            planes.append(arr)
        frames.append(planes)
    def fill(n, f):
        dst = f.copy()
        for p, samples in enumerate(frames[n]): np.copyto(np.asarray(dst[p]), samples)
        dst.props['Fixture'] = b'neo-fft-p1-v1\x00payload'
        dst.props['Sequence'] = n
        dst.props['_ColorRange'] = 1
        dst.props['_ChromaLocation'] = 0
        return dst
    node = vs.core.std.ModifyFrame(base, clips=base, selector=fill)
    digest = hashlib.sha256(b''.join(a.tobytes() for f in frames for a in f)).hexdigest()
    return node, frames, digest

def cases():
    out = []
    def add(algorithm, group, params, **kw):
        out.append(dict(id=f'{algorithm}-{group}-{len(out):04}', algorithm=algorithm,
                        format=kw.pop('format', 'gray'), bits=kw.pop('bits', 32), params=params, **kw))
    for a in ('FFT3D', 'DFTTest'):
        for family in ('gray', '444', '422', '420', 'rgb'):
            for bits in (8, 10, 12, 14, 16, 32):
                add(a, 'baseline', {}, format=family, bits=bits, frames=7)
                add(a, 'plane-copy', {'planes': [0]}, format=family, bits=bits)
    for i in range(18):
        bw, bh = ((8,8), (8,6), (9,7), (32,32))[i%4]
        ow, oh = ((0,0), (-7,-3), (bw//2,bh//2), (bw//3,1))[i%4]
        add('FFT3D','spatial',dict(bw=bw,bh=bh,ow=ow,oh=oh,wintype=i%3,sigma=(0,2,8)[i%3],
                                 beta=1+i%2,degrid=(0,.5,1)[(i//3)%3]), width=130,height=98)
    for win in range(3):
        for degrid in (0,1):
            add('FFT3D','windows',dict(bw=9,bh=7,ow=4,oh=2,wintype=win,degrid=degrid,sigma=8))
            add('FFT3D','windows-yuv',dict(bw=8,bh=6,ow=4,oh=2,wintype=win,degrid=degrid),format='420',bits=16)
    for mode in (0,1):
        for ftype in range(5):
            for zmean in (False,True):
                add('DFTTest','filter',dict(sbsize=9 if mode==0 else 8,smode=mode,sosize=6,ftype=ftype,zmean=zmean,
                    sigma=32 if ftype<2 else .25,sigma2=1.5,pmin=2,pmax=100))
        for win in range(12):
            add('DFTTest','window',dict(sbsize=9 if mode==0 else 8,smode=mode,sosize=4,swin=win,
                ftype=win%5,zmean=bool(win%2),sigma=8 if win%5<2 else .5))
            add('DFTTest','temporal-window',dict(sbsize=9 if mode==0 else 8,smode=mode,sosize=4,twin=win))
            add('DFTTest','window-yuv',dict(sbsize=9 if mode==0 else 8,smode=mode,sosize=4,swin=win,
                zmean=False,sigma=8),format='444',bits=16)
    for mode in (0,1):
        for exponent in (.5,1,2,.49995,.50005,.99995,1.00005):
            for beta in (np.nextafter(np.float32(exponent),np.float32(0)),np.nextafter(np.float32(exponent),np.float32(3))):
                add('DFTTest','beta',dict(sbsize=3 if mode==0 else 8,smode=mode,sosize=4,f0beta=float(beta),zmean=False))
    for b, o in ((8,0),(8,4),(8,6),(9,0),(9,4),(16,0),(16,8),(16,12)):
        add('DFTTest','grid',dict(sbsize=b,sosize=o),width=130,height=98)
    for b in (1,3,9): add('DFTTest','center',dict(sbsize=b,smode=0,ftype=2,sigma=.5,zmean=False))
    for a in ('FFT3D','DFTTest'):
        for planes in ([],[1],[0,1,2],[2,0,2,1]): add(a,'planes',dict(planes=planes),format='420',bits=32)
    for ftype in (0,1,2):
        for sigma in ((0,8,32) if ftype<2 else (0,.5,1)):
            add('DFTTest','values',dict(ftype=ftype,sigma=sigma,zmean=False,sbsize=8,sosize=4))
    for mode in (0,1):
        add('DFTTest','kaiser-zero',dict(sbsize=9 if mode==0 else 8,smode=mode,sosize=4,swin=4,twin=4,sbeta=0,tbeta=0))
    return out

def supplemental_cases():
    """Explicit spec coverage additions; retain v1 calibration inputs unchanged."""
    out = []
    def add(a, group, params, **kw):
        out.append(dict(id=f'{a}-supplement-{group}-{len(out):04}', algorithm=a,
                        format=kw.pop('format', 'gray'), bits=kw.pop('bits', 32), params=params, **kw))
    for a in ('FFT3D', 'DFTTest'):
        for family in ('422', '420'):
            for bits in (16, 32):
                add(a, 'nondivisible', {}, format=family, bits=bits, width=260, height=196, frames=7)
        add(a, 'signed-extremes', {}, pattern=7)
        for pattern in (8, 9, 10, 11): add(a, 'directional', {}, pattern=pattern)
    for mode in (0, 1):
        for beta in (.5, 1., 2.):
            add('DFTTest', 'exact-beta', dict(sbsize=3 if mode==0 else 8, smode=mode, sosize=4,
                                           f0beta=beta, zmean=False))
    add('DFTTest', 'zero-dc-no-mean', dict(sbsize=4, sosize=0, swin=6, zmean=False))
    for ftype in (3, 4):
        for pmin, pmax in ((0., 1.), (2., 10.)):
            add('DFTTest', 'bands', dict(sbsize=8, sosize=4, ftype=ftype, sigma=.25,
                                       sigma2=1.5, pmin=pmin, pmax=pmax, zmean=False))
    # First safe reflection sizes for B=8/O=0, and B=1 overlap-add.
    add('FFT3D', 'first-reflection', dict(bw=8, bh=8, ow=0, oh=0), width=13, height=13)
    add('DFTTest', 'first-reflection', dict(sbsize=1, sosize=0), width=2, height=2)
    return out
