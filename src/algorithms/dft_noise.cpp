#include "algorithms/dft_noise.hpp"
#include "algorithms/plan.hpp"

namespace neo_fft {
DFTNoise::DFTNoise(const DFTConfig& c)
    : locations(c.locations), temporal_size(c.tbsize), block_size(c.block),
      fft_(c.tbsize,c.block,c.block), zmean_(c.zmean) {
  validate(c);
  require(c.ftype < 2 && !locations.empty(), "inactive DFTTest sample model");
  auto sample = dft_window_3d(c.tbsize,c.block,0,0,c.swin,c.twin,c.sbeta,c.tbeta);
  const auto output = dft_window_3d(c.tbsize,c.block,c.mode == 0 ? 0 : c.overlap,c.mode,c.swin,c.twin,c.sbeta,c.tbeta);
  calibration_ = finite(finite((1.0f / float(locations.size())) * finite(sample.wscale/output.wscale)) * c.alpha.value_or(c.ftype == 0 ? 5.0f : 7.0f));
  window_ = std::move(sample.h);
  auto input = window_;
  for (float& v : input) v = finite(255.0f*v);
  grid_ = buffer<std::complex<float>>(fft_.bins());
  fft_.forward(input.data(),grid_.data());
  for (auto v : grid_) { finite(v.real()); finite(v.imag()); }
  require(!zmean_ || grid_[0].real()!=0, "DFTTest sample template DC must be nonzero");
}
void DFTNoise::prepare(const Gather& gather) const {
  if (model_.get()) return;
  auto input = buffer<float>(fft_.samples());
  auto spectrum = buffer<std::complex<float>>(fft_.bins());
  auto sum = buffer<float>(fft_.bins());
  const auto slice = mul_size(std::size_t(block_size),std::size_t(block_size));
  for (const auto& location : locations) {
    for (int z=0;z<temporal_size;++z)
      gather(location,z,{input.data()+std::size_t(z)*slice,slice});
    for (std::size_t k=0;k<input.size();++k) input[k]=finite(finite(input[k])*window_[k]);
    fft_.forward(input.data(),spectrum.data());
    const float ratio=zmean_ ? finite(spectrum[0].real()/grid_[0].real()) : 0;
    for (std::size_t k=0;k<sum.size();++k) {
      const float re=finite(spectrum[k].real()-(zmean_ ? finite(ratio*grid_[k].real()) : 0));
      const float im=finite(spectrum[k].imag()-(zmean_ ? finite(ratio*grid_[k].imag()) : 0));
      sum[k]=finite(sum[k]+finite(finite(re*re)+finite(im*im)));
    }
  }
  for (float& v : sum) v=finite(v*calibration_);
  model_.publish(std::move(sum));
}
} // namespace neo_fft
