#pragma once
#include "base/checked.hpp"
#include <array>
namespace neo_fft {
inline void validate_temporal(int size,int mode,int overlap) {
  require(mode==0 || mode==1,"DFTTest tmode outside 0..1");
  require(size>=1 && size<=15,"DFTTest tbsize outside 1..15");
  require(mode!=0 || size%2==1,"DFTTest centered tbsize must be odd");
  if(mode==1) {
    require(overlap>=0 && overlap<size,"DFTTest invalid tosize");
    require(overlap<=size/2 || size%(size-overlap)==0,"DFTTest temporal heavy overlap requires divisible step");
  }
}
struct TemporalBlock {
  std::int64_t start;
  int target;
  std::array<int,15> slots{};
};
inline std::int64_t floor_div(std::int64_t value,int divisor) {
  require(divisor>0,"invalid temporal divisor");
  return value/divisor-(value%divisor<0);
}
// Absolute lattice, bounded by 15 covering blocks even at INT_MAX indices.
inline std::vector<TemporalBlock> temporal_blocks(int frame,int length,int size,int overlap) {
  validate_temporal(size,1,overlap);
  require(length>=size && frame>=0 && frame<length,"invalid temporal clip/index");
  const int step=size-overlap;
  const std::int64_t anchor=-std::max(step,overlap),n=frame;
  const auto first=std::max(std::int64_t(0),floor_div(n-size-anchor,step)+1);
  const auto last=floor_div(n-anchor,step);
  std::vector<TemporalBlock> blocks;
  blocks.reserve(std::size_t(last-first+1));
  for(auto j=first;j<=last;++j) {
    const auto start=anchor+j*step;
    TemporalBlock b{start,int(n-start),{}};
    for(int z=0;z<size;++z)b.slots[z]=int(std::clamp(start+z,std::int64_t(0),std::int64_t(length)-1));
    blocks.push_back(b);
  }
  return blocks;
}
} // namespace neo_fft
