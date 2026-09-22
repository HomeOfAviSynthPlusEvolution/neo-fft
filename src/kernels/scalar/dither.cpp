#include "kernels/dither.hpp"
namespace neo_fft {
namespace {
std::uint32_t mix(std::uint32_t v) {
  v ^= v>>16; v *= 0x7feb352du; v ^= v>>15; v *= 0x846ca68bu; return v^(v>>16);
}
}
std::uint32_t dither_hash(std::uint32_t seed,std::uint32_t n,std::uint32_t p,std::uint32_t y,std::uint32_t x) {
  return mix(mix(mix(mix(mix(seed^0xa511e9b3u)^n)^p)^y)^x);
}
void dither_scalar(span2d::Plane<const float> source, span2d::Plane<std::uint8_t> dst,
                   int mode, int seed, int frame, int plane,DitherNoise noise) {
  require(mode>=1 && seed>=0 && frame>=0 && plane>=0,"invalid dither coordinates or mode");
  require(source.width()==dst.width() && source.height()==dst.height(),"dither dimensions differ");
  auto current=buffer<float>(source.width()), next=buffer<float>(source.width());
  const float scale=float(mode-1)+.5f, off=scale*.5f;
  alignas(64) float random[64];
  for (int y=0;y<source.height();++y) {
    std::fill(next.begin(),next.end(),0);
    for (int x=0;x<source.width();++x) {
      const float e=finite(source.row_ptr(y)[x]);
      float v;
      if (mode==1) v=finite(finite(e+current[x])+.5f);
      else {
        if(x%64==0)noise(random,std::size_t(std::min(64,source.width()-x)),seed,frame,plane,y,x);
        const float u=random[x%64];
        v=finite(finite(finite(finite(e+u*scale)-off)+current[x])+.5f);
      }
      const auto d=static_cast<std::uint8_t>(v<=0 ? 0 : v>=255 ? 255 : v);
      dst.row_ptr(y)[x]=d;
      const float error=finite(e-float(d));
      if (x>0) next[x-1]=finite(next[x-1]+error*.1875f);
      next[x]=finite(next[x]+error*.3125f);
      if (x+1<source.width()) {
        current[x+1]=finite(current[x+1]+error*.4375f);
        next[x+1]=finite(next[x+1]+error*.0625f);
      }
    }
    current.swap(next);
  }
}
} // namespace neo_fft
