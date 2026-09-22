#include "kernels/rows.hpp"
#include "kernels/dither.hpp"
#include <cstring>
namespace neo_fft {
void copy_row_scalar(const void* source,void* dst,std::size_t bytes) {if(bytes)std::memcpy(dst,source,bytes);}
void dither_noise_scalar(float* dst,std::size_t count,std::uint32_t seed,std::uint32_t frame,std::uint32_t plane,std::uint32_t y,std::uint32_t x) {
  for(std::size_t i=0;i<count;++i)dst[i]=float(dither_hash(seed,frame,plane,y,x+std::uint32_t(i))>>8)*0x1p-24f;
}
#if !NEO_FFT_ENABLE_HIGHWAY
CopyRow select_copy_row(int) {return copy_row_scalar;}
DitherNoise select_dither_noise(int) {return dither_noise_scalar;}
#endif
} // namespace neo_fft
