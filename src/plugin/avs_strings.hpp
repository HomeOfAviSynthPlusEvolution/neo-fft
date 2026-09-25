#pragma once
#include <charconv>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace neo_fft::plugin::avs {

inline bool accepts_dft_text_array(std::string_view name) noexcept {
  return name == "nlocation" || name == "slocation" || name == "ssx" || name == "ssy" || name == "sst";
}

inline bool decimal_token(std::string_view token, bool floating) noexcept {
  std::size_t i = 0;
  const auto sign = [&] {
    if (i < token.size() && (token[i] == '+' || token[i] == '-')) ++i;
  };
  const auto digits = [&] {
    const auto start = i;
    while (i < token.size() && token[i] >= '0' && token[i] <= '9') ++i;
    return i - start;
  };
  sign();
  auto count = digits();
  if (floating && i < token.size() && token[i] == '.') {
    ++i;
    count += digits();
  }
  if (!count) return false;
  if (floating && i < token.size() && (token[i] == 'e' || token[i] == 'E')) {
    ++i;
    sign();
    if (!digits()) return false;
  }
  return i == token.size();
}

// Takes a validated decimal_token; an unavailable floating overload is never compiled.
bool parse_double_token(std::string_view token, double& value);

template<class T>
T parse_number_token(const char* name, std::string_view token) {
  static_assert(std::is_same_v<T, int> || std::is_same_v<T, double>);
  const auto invalid = [&] {
    return std::invalid_argument(std::string("DFTTest ") + name + ": invalid decimal string token");
  };
  if (!decimal_token(token, std::is_floating_point_v<T>)) throw invalid();
  T value{};
  if constexpr (std::is_same_v<T, int>) {
    // from_chars does not accept a leading '+'. Syntax was checked above.
    if (token.front() == '+') token.remove_prefix(1);
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    if (result.ec != std::errc{} || result.ptr != token.data() + token.size()) throw invalid();
  } else {
    if (!parse_double_token(token, value)) throw invalid();
  }
  return value;
}

// Creation-time adapter only: preserve the common numeric-array validation.
template <class T>
std::vector<T> parse_number_list(const char* name, std::string_view text) {
  // Upstream AVS lists allow whitespace, commas and colons, also mixed.
  constexpr std::string_view separators = " \t\r\n\v\f,:";
  std::vector<T> values;
  auto first = text.find_first_not_of(separators);
  while (first != std::string_view::npos) {
    const auto end = text.find_first_of(separators, first);
    values.push_back(parse_number_token<T>(name, text.substr(first, end == std::string_view::npos ? end : end - first)));
    if (end == std::string_view::npos) break;
    first = text.find_first_not_of(separators, end);
  }
  return values;
}

} // namespace neo_fft::plugin::avs
