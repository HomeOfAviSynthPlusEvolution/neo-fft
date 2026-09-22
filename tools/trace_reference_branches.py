"""Windows x64 black-box branch evidence via one-shot child-process breakpoints.

Uses only PE/COFF symbols, never reference source. The DLL on disk is untouched.
Each breakpoint is restored at its first hit before continuing the original code.
Only an explicitly spawned single-worker probe is debugged, never user processes.
"""
import argparse
import ctypes as C
from ctypes import wintypes as W
import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess
import sys
import time

def child(a):
 sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tests/vapoursynth'))
 from fixtures import environment
 import vapoursynth as vs
 environment(vs);core=vs.core;core.num_threads=1
 core.std.LoadPlugin(path=str(a.plugin.resolve()))
 src=core.std.BlankClip(width=64,height=48,length=7,format=vs.GRAYS,color=[.25])
 if a.algorithm=='FFT3D':
  kw=dict(bw=8,bh=8,bt=a.temporal,sharpen=.3,dehalo=.2,fft_backend='pocketfft',ncpu=1)
  if a.mode=='analytic':kw.update(sigma2=3)
  if a.mode=='sample':kw.update(pfactor=1,px=1,py=1,pframe=0)
  core.neo_fft3d.FFT3D(src,opt=a.opt,**kw).get_frame(3)
 else:
  kw=dict(sbsize=5,sosize=2,tbsize=a.temporal,fft_backend='pocketfft',threads=1,fft_threads=1)
  if a.mode=='sample':kw.update(nlocation=[0,0,0,0])
  elif a.mode=='analytic':kw.update(slocation=[0,2,1,12])
  core.neo_dfttest.DFTTest(src,opt=a.opt,**kw).get_frame(3)

def trace(a):
 if sys.platform!='win32' or C.sizeof(C.c_void_p)!=8:raise RuntimeError('Windows x64 required')
 data=a.plugin.read_bytes();pe=struct.unpack_from('<I',data,0x3c)[0]
 if data[pe:pe+4]!=b'PE\0\0' or struct.unpack_from('<H',data,pe+4)[0]!=0x8664:raise RuntimeError('x64 PE required')
 preferred=struct.unpack_from('<Q',data,pe+48)[0];size=struct.unpack_from('<I',data,pe+80)[0]
 symbols={}
 for line in subprocess.check_output([str(a.nm),'-C',str(a.plugin)],text=True).splitlines():
  fields=line.split(maxsplit=2)
  if len(fields)!=3 or fields[1]!='T':continue
  name=fields[2]
  if a.algorithm=='FFT3D':selected=re.search(r'\b(?:Apply[23]D|Sharpen)(?:[2345])?_(?:C|Hwy)',name)
  else:selected=('filter_scalar<' in name or re.search(r'process_(?:spatial|temporal)_scalar<',name) or re.search(r'::Process(?:Spatial|Temporal)Highway<',name)) and 'Pipeline' not in name
  if selected:
   rva=int(fields[0],16)-preferred
   if not 0<=rva<size:raise RuntimeError('symbol outside image')
   symbols[rva]=name
 if not symbols:raise RuntimeError('no traceable symbols')
 k=C.WinDLL('kernel32',use_last_error=True)
 def fn(name,args,result=W.BOOL):
  f=getattr(k,name);f.argtypes=args;f.restype=result;return f
 close=fn('CloseHandle',[W.HANDLE]);read=fn('ReadProcessMemory',[W.HANDLE,C.c_void_p,C.c_void_p,C.c_size_t,C.c_void_p])
 write=fn('WriteProcessMemory',[W.HANDLE,C.c_void_p,C.c_void_p,C.c_size_t,C.c_void_p])
 flush=fn('FlushInstructionCache',[W.HANDLE,C.c_void_p,C.c_size_t])
 wait=fn('WaitForDebugEvent',[C.c_void_p,W.DWORD]);cont=fn('ContinueDebugEvent',[W.DWORD,W.DWORD,W.DWORD])
 open_thread=fn('OpenThread',[W.DWORD,W.BOOL,W.DWORD],W.HANDLE)
 getctx=fn('GetThreadContext',[W.HANDLE,C.c_void_p]);setctx=fn('SetThreadContext',[W.HANDLE,C.c_void_p])
 pathfn=fn('GetFinalPathNameByHandleW',[W.HANDLE,W.LPWSTR,W.DWORD,W.DWORD],W.DWORD)
 create=fn('CreateProcessW',[W.LPCWSTR,W.LPWSTR,C.c_void_p,C.c_void_p,W.BOOL,W.DWORD,C.c_void_p,W.LPCWSTR,C.c_void_p,C.c_void_p])
 terminate=fn('TerminateProcess',[W.HANDLE,W.UINT])
 def check(ok):
  if not ok:raise C.WinError(C.get_last_error())
 command=subprocess.list2cmdline([sys.executable,str(Path(__file__).resolve()),'--child','--plugin',str(a.plugin.resolve()),'--algorithm',a.algorithm,'--mode',a.mode,'--temporal',str(a.temporal),'--opt',str(a.opt)])
 startup=C.create_string_buffer(104);struct.pack_into('<I',startup,0,104);process=C.create_string_buffer(24)
 check(create(None,C.create_unicode_buffer(command),None,None,False,2|0x08000000,None,None,startup,process))
 handle,main_thread,pid,tid=struct.unpack_from('<QQII',process);close(main_thread)
 event=C.create_string_buffer(176);armed={};hits=[];installed=False;exited=False;deadline=time.monotonic()+90
 try:
  while time.monotonic()<deadline:
   if not wait(event,1000):
    if C.get_last_error()==121:continue
    check(False)
   code,pid,tid=struct.unpack_from('<III',event);status=0x10002
   if code==6:
    file,base=struct.unpack_from('<QQ',event,16)
    if file:
     path=C.create_unicode_buffer(32768);length=pathfn(file,path,len(path),0)
     if length and Path(path.value).name.lower()==a.plugin.name.lower():
      for rva,name in symbols.items():
       address=base+rva;original=C.create_string_buffer(1);check(read(handle,address,original,1,None))
       check(write(handle,address,C.c_char_p(b'\xcc'),1,None));check(flush(handle,address,1))
       armed[address]=(original.raw,name)
      installed=True
     close(file)
   elif code==3:
    for offset in (16,24,32):
     event_handle=struct.unpack_from('<Q',event,offset)[0]
     if event_handle:close(event_handle)
   elif code==2:
    event_handle=struct.unpack_from('<Q',event,16)[0]
    if event_handle:close(event_handle)
   elif code==1:
    exception=struct.unpack_from('<I',event,16)[0];address=struct.unpack_from('<Q',event,32)[0]
    if exception==0x80000003 and address in armed:
     original,name=armed.pop(address);check(write(handle,address,C.c_char_p(original),1,None));check(flush(handle,address,1))
     thread=open_thread(0x1a,False,tid);check(thread)
     try:
      storage=C.create_string_buffer(1232+16);aligned=(C.addressof(storage)+15)&~15
      C.c_uint32.from_address(aligned+48).value=0x100001;check(getctx(thread,aligned))
      C.c_uint64.from_address(aligned+248).value=address;check(setctx(thread,aligned))
     finally:close(thread)
     hits.append(name)
    elif exception!=0x80000003:status=0x80010001
   elif code==5:
    exit_code=struct.unpack_from('<I',event,16)[0];exited=True
   check(cont(pid,tid,status))
   if exited:break
  if not exited:raise RuntimeError('probe timed out')
  if exit_code or not installed or not hits:raise RuntimeError(f'probe failed: exit={exit_code:#x}, installed={installed}, hits={hits}')
 finally:
  if not exited:terminate(handle,1)
  close(handle)
 report=dict(algorithm=a.algorithm,mode=a.mode,temporal=a.temporal,opt=a.opt,sha256=hashlib.sha256(data).hexdigest(),executed_symbols=hits)
 a.output.write_text(json.dumps(report,indent=2));print(json.dumps(report))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--plugin',type=Path,required=True);p.add_argument('--algorithm',choices=['FFT3D','DFTTest'],required=True);p.add_argument('--mode',choices=['uniform','analytic','sample'],default='analytic');p.add_argument('--temporal',type=int,default=1);p.add_argument('--opt',type=int,default=0);p.add_argument('--child',action='store_true');p.add_argument('--nm',type=Path);p.add_argument('--output',type=Path);args=p.parse_args()
 if args.child:child(args)
 else:trace(args)
