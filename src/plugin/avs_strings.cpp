#include "avs_strings.hpp"

#ifndef NEO_FFT_HAS_DOUBLE_FROM_CHARS
#include <cerrno>
#include <cstdlib>
#include <locale.h>
#include <limits>
#include <new>
#if defined(__APPLE__) || defined(__FreeBSD__)
#include <xlocale.h>
#endif
#endif

namespace neo_fft::plugin::avs {
#ifndef NEO_FFT_HAS_DOUBLE_FROM_CHARS
namespace {
struct NumericLocale {
#ifdef _WIN32
  _locale_t handle = _create_locale(LC_NUMERIC, "C");
  ~NumericLocale() { _free_locale(handle); }
#else
  locale_t handle = newlocale(LC_NUMERIC_MASK, "C", nullptr);
  ~NumericLocale() { freelocale(handle); }
#endif
  NumericLocale() {
    if (!handle) throw std::bad_alloc();
  }
  NumericLocale(const NumericLocale&) = delete;
  NumericLocale& operator=(const NumericLocale&) = delete;
};

struct PreserveErrno {
  int saved = errno;
  ~PreserveErrno() { errno = saved; }
};
} // namespace
#endif

bool parse_double_token(std::string_view token, double& value) {
#ifdef NEO_FFT_HAS_DOUBLE_FROM_CHARS
  if (token.front() == '+') token.remove_prefix(1);
  const auto result = std::from_chars(token.data(), token.data() + token.size(), value,
                                    std::chars_format::general);
  return result.ec == std::errc{} && result.ptr == token.data() + token.size() && std::isfinite(value);
#else
  const PreserveErrno preserve_errno;
  static const NumericLocale locale;
  const std::string text(token);
  char* end = nullptr;
  errno = 0;
#ifdef _WIN32
  const double parsed = _strtod_l(text.c_str(), &end, locale.handle);
#else
  const double parsed = strtod_l(text.c_str(), &end, locale.handle);
#endif
  // Accept representable underflow (including rounding up to the minimum normal).
  // Directed rounding can make overflow return finite DBL_MAX instead of infinity.
  const bool range_error = errno == ERANGE &&
      (parsed == 0 || std::abs(parsed) > std::numeric_limits<double>::min());
  if (end != text.c_str() + text.size() || !std::isfinite(parsed) || range_error)
    return false;
  value = parsed;
  return true;
#endif
}
} // namespace neo_fft::plugin::avs
