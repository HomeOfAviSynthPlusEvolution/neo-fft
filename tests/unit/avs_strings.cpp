#include "plugin/avs_strings.hpp"
#include "../test.hpp"
#include <limits>

using namespace neo_fft::plugin::avs;

struct CommaDecimal : std::numpunct<char> {
  char do_decimal_point() const override { return ','; }
};

template<class T>
void rejects_token(const char* name, const char* text) {
  try {
    parse_number_list<T>(name, text);
  } catch (const std::invalid_argument&) {
    return;
  }
  throw std::runtime_error(std::string("Unexpectedly accepted numeric string: ") + text);
}

int main() {
  try {
    CHECK(parse_number_list<int>("nlocation", " \t\r\n\v\f").empty());
    CHECK(parse_number_list<double>("slocation", "").empty());
    CHECK((parse_number_list<int>("nlocation", " +1\t-2\r\n0 03 ") == std::vector<int>{1, -2, 0, 3}));
    CHECK((parse_number_list<int>("nlocation", "-2147483648 2147483647") ==
           std::vector<int>{std::numeric_limits<int>::min(), std::numeric_limits<int>::max()}));
    CHECK((parse_number_list<double>("ssx", "0 +.5 1. 2.5e+1 1E-2 -0") ==
           std::vector<double>{0, .5, 1, 25, .01, 0}));
    for (const char* text : {"1.5", "1e2", "1+2", "0x10", "2147483648", "-2147483649",
                             "1,2", "[1 2]", "1 bad", "+", "--1", "nan", "inf"})
      rejects_token<int>("nlocation", text);
    for (const char* text : {".", "+.", "1e+", "1e-", "1+2", "0x1p2", "1e", "1e999", "1,5", "[0 1]", "0:1",
                             "nan", "NaN", "inf", "-infinity", "1 2junk"})
      rejects_token<double>("slocation", text);
    const auto saved = std::locale();
    std::locale::global(std::locale(saved, new CommaDecimal));
    const auto values = parse_number_list<double>("ssy", "0 .5 1 2.25");
    std::locale::global(saved);
    CHECK((values == std::vector<double>{0, .5, 1, 2.25}));
    try {
      parse_number_list<int>("nlocation", "1.5");
      CHECK(false);
    } catch (const std::invalid_argument& error) {
      CHECK(std::string(error.what()).find("nlocation") != std::string::npos);
    }
    for (const char* name : {"nlocation", "slocation", "ssx", "ssy", "sst"})
      CHECK(accepts_dft_text_array(name));
    CHECK(!accepts_dft_text_array("planes"));
    // Exercise the older-libc++ conversion fallback even on libraries that
    // provide floating from_chars; both paths keep the same lexical contract.
    for (const char* token : {"0", "-0", "+.5", "1.", "2.5e+1", "1E-2", "1e-300", "1e300",
                              "4.9406564584124654e-324", "2.2250738585072014e-308"})
      CHECK((parse_number_token<double, false>("sst", token) == parse_number_token<double>("sst", token)));
    for (const char* token : {"0x1p2", "1e+", ".", "nan", "inf", "1e999", "1e-999", "1,5"}) {
      rejects([&] { parse_number_token<double>("sst", token); });
      rejects([&] { parse_number_token<double, false>("sst", token); });
    }
    std::cout << "AVS numeric strings: strict tokens, bounds and locale passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
