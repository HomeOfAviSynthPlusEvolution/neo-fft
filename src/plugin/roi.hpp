#pragma once
#include "algorithms/plan.hpp"
#include <cstring>

namespace neo_fft::plugin {
struct ROI {
  int left = 0, top = 0, width = 0, height = 0;
  bool interlaced = false;
  int row(int y) const noexcept {
    return interlaced ? (y < height / 2 ? 2*y : 2*(height-1-y)+1) : y;
  }
};
inline ROI make_roi(int width, int height, int sub_w, int sub_h, const FFT3DConfig& c) {
  const int ax = 1 << sub_w, ay = 1 << sub_h;
  require(c.left % ax == 0 && c.right % ax == 0 && c.top % ay == 0 && c.bottom % ay == 0,
          "selected chroma ROI margins must align to subsampling");
  const int l = c.left / ax, r = c.right / ax, t = c.top / ay, b = c.bottom / ay;
  require(l < width && r < width-l && t < height && b < height-t, "ROI must be nonempty");
  ROI roi{l,t,width-l-r,height-t-b,c.interlaced};
  require(!c.interlaced || roi.height % 2 == 0, "interlaced ROI height must be even");
  return roi;
}
// A request owns packed pixels; all reflection stays inside this working image.
template<class T> struct PackedROI {
  std::vector<T> pixels;
  span2d::Plane<T> view;
  PackedROI(span2d::Plane<const T> source, const ROI& roi)
      : pixels(mul_size(std::size_t(roi.width),std::size_t(roi.height))) {
    checked_subplane(source,roi.left,roi.top,roi.width,roi.height);
    const auto pitch = mul_size(std::size_t(roi.width),sizeof(T));
    view = checked_plane(pixels.data(),roi.width,roi.height,pitch,mul_size(pixels.size(),sizeof(T)));
    for (int y=0;y<roi.height;++y)
      std::memcpy(view.row_ptr(y),source.row_ptr(roi.top+roi.row(y))+roi.left,pitch);
  }
  void write(span2d::Plane<T> dst, const ROI& roi) const {
    for (int y=0;y<roi.height;++y)
      std::memcpy(dst.row_ptr(roi.top+roi.row(y))+roi.left,view.row_ptr(y),std::size_t(roi.width)*sizeof(T));
  }
};
} // namespace neo_fft::plugin
