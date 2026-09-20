#pragma once
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#define CHECK(x)                                                                                                       \
  do {                                                                                                                 \
    if (!(x))                                                                                                          \
      throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " #x);                      \
  } while (0)
inline void check_near(double a, double b, double tolerance = 2e-5) {
  CHECK(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tolerance);
}
template <class F>
void rejects(F f) {
  bool caught = false;
  try {
    f();
  } catch (const std::exception&) {
    caught = true;
  }
  CHECK(caught);
}
