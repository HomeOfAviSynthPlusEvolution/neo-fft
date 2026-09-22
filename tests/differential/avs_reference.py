"""Small AviSynth C API pixel reader; uses public host API only."""
import ctypes as C
import os
import numpy as np
class Data(C.Union):
    _fields_=[('pointer',C.c_void_p),('string',C.c_char_p),('integer',C.c_int),('floating',C.c_float)]
class Value(C.Structure):
    _fields_=[('type',C.c_short),('size',C.c_short),('data',Data)]
class AviSynth:
    def __init__(self,dll,fftw):
        self.search=os.add_dll_directory(str(fftw.parent))
        self.fftw=C.WinDLL(str(fftw))
        self.dll=C.WinDLL(str(dll))
        signatures={
          'avs_create_script_environment':(C.c_void_p,[C.c_int]),
          'avs_delete_script_environment':(None,[C.c_void_p]),
          'avs_invoke':(Value,[C.c_void_p,C.c_char_p,Value,C.c_void_p]),
          'avs_take_clip':(C.c_void_p,[Value,C.c_void_p]),
          'avs_release_value':(None,[Value]),'avs_release_clip':(None,[C.c_void_p]),
          'avs_get_frame':(C.c_void_p,[C.c_void_p,C.c_int]),
          'avs_clip_get_error':(C.c_char_p,[C.c_void_p]),
          'avs_release_video_frame':(None,[C.c_void_p]),
          'avs_get_read_ptr_p':(C.c_void_p,[C.c_void_p,C.c_int]),
          'avs_get_pitch_p':(C.c_int,[C.c_void_p,C.c_int]),
          'avs_get_row_size_p':(C.c_int,[C.c_void_p,C.c_int]),
          'avs_get_height_p':(C.c_int,[C.c_void_p,C.c_int])}
        for name,(result,args) in signatures.items():
            fn=getattr(self.dll,name);fn.restype=result;fn.argtypes=args
        self.env=self.dll.avs_create_script_environment(8)
        if not self.env:raise RuntimeError('AviSynth environment unavailable')
        self.clips=[]
    def evaluate(self,script):
        raw=script.encode();arg=Value(ord('s'),0,Data(string=raw))
        value=self.dll.avs_invoke(self.env,b'Eval',arg,None)
        try:
            if value.type==ord('e'):raise RuntimeError(value.data.string.decode(errors='replace'))
            if value.type!=ord('c'):raise RuntimeError('script did not return a clip')
            clip=self.dll.avs_take_clip(value,self.env);self.clips.append(clip);return clip
        finally:self.dll.avs_release_value(value)
    def frame(self,clip,n,bits,family):
        frame=self.dll.avs_get_frame(clip,n)
        if not frame:raise RuntimeError(self.dll.avs_clip_get_error(clip))
        try:
            planes=[1] if family=='gray' else ([32,64,128] if family=='rgb' else [1,2,4])
            dtype=np.float32 if bits==32 else (np.uint8 if bits==8 else np.uint16)
            result=[]
            for plane in planes:
                ptr=self.dll.avs_get_read_ptr_p(frame,plane);pitch=self.dll.avs_get_pitch_p(frame,plane)
                size=self.dll.avs_get_row_size_p(frame,plane);height=self.dll.avs_get_height_p(frame,plane)
                raw=b''.join(C.string_at(ptr+y*pitch,size) for y in range(height))
                result.append(np.frombuffer(raw,dtype=dtype).reshape(height,-1).copy())
            return result
        finally:self.dll.avs_release_video_frame(frame)
    def drop(self,clip):
        self.dll.avs_release_clip(clip);self.clips.remove(clip)
    def close(self):
        for clip in reversed(self.clips):self.dll.avs_release_clip(clip)
        self.clips=[]
        if self.env:self.dll.avs_delete_script_environment(self.env);self.env=None
        self.search.close()
