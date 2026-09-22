#pragma once
#include <cstddef>
namespace neo_fft {
// Creation-time table rows. Coordinates/knots are validated by the builders.
// Exact libm exp/pow and ordered reductions stay canonical; bulk arithmetic is
// selected with the same opt policy as filtering. No padding is required.
// maximum consumes finite nonnegative tables. window carries the ordered
// binary32 energy sum between rows after binary64 tensor multiplication.
struct TableKernels {
  void (*radius)(float*,const float*,std::size_t,float squared_other,float dimensions);
  void (*analytic)(float*,const float*,std::size_t,const float* sigmas,float norm);
  void (*curve)(float*,const float*,std::size_t,const float*,const float*,std::size_t,float divisor);
  void (*product)(float*,const float*,std::size_t,float factor,float divisor);
  void (*weights)(float*,const float*,std::size_t,float squared_other,float cutoff_squared);
  void (*enhancement)(float*,float*,const float*,std::size_t,float squared_other,float cutoff,float hr);
  void (*divide)(float*,std::size_t,float divisor);
  float (*window)(float*,const double*,std::size_t,double factor,double normalization,float energy);
  float (*maximum)(const float*,std::size_t);
  void (*check_scale)(const float*,std::size_t,float);
};
const TableKernels& table_scalar();
TableKernels select_table(int opt);
}
