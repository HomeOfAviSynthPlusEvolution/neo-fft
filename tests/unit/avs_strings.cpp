#include "plugin/avs_strings.hpp"
#include "../test.hpp"
#include <cerrno>
#include <cfenv>
#include <clocale>
#include <limits>
#include <locale>
#include <utility>

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
    CHECK(parse_number_list<int>("nlocation", " ,:: \t, \r\n").empty());
    CHECK(parse_number_list<double>("slocation", ":,,").empty());
    for (const char* separator : {" ", ",", ":", ",: \t\r\n\v\f:,"}) {
      const std::string sep=separator;
      CHECK((parse_number_list<int>("nlocation", sep+"+0"+sep+"1"+sep+"2"+sep+"3"+sep) ==
             std::vector<int>{0,1,2,3}));
      CHECK((parse_number_list<double>("slocation", sep+"0"+sep+"+.5e-1"+sep+"1"+sep+"2.5E+1"+sep) ==
             std::vector<double>{0,.05,1,25}));
    }
    // Commas delimit values; they never act as locale-specific decimal points.
    CHECK((parse_number_list<double>("ssx", "1,5") == std::vector<double>{1,5}));
    CHECK((parse_number_list<int>("nlocation", " +1\t-2\r\n0 03 ") == std::vector<int>{1, -2, 0, 3}));
    CHECK((parse_number_list<int>("nlocation", "-2147483648 2147483647") ==
           std::vector<int>{std::numeric_limits<int>::min(), std::numeric_limits<int>::max()}));
    CHECK((parse_number_list<double>("ssx", "0 +.5 1. 2.5e+1 1E-2 -0") ==
           std::vector<double>{0, .5, 1, 25, .01, 0}));
    for (const char* text : {"1.5", "1e2", "1+2", "0x10", "2147483648", "-2147483649",
                             "1;2", "1,:bad", "1e,2", "[1 2]", "1 bad", "+", "--1", "nan", "inf"})
      rejects_token<int>("nlocation", text);
    for (const char* text : {".", "+.", "1e+", "1e-", "1+2", "0x1p2", "1e", "1e999", "1;5", "1e+:2", "0,:junk", "[0 1]",
                             "nan", "NaN", "inf", "-infinity", "1 2junk"})
      rejects_token<double>("slocation", text);
    const auto saved = std::locale();
    std::locale::global(std::locale(saved, new CommaDecimal));
    const auto values = parse_number_list<double>("ssy", "0:.5,1 2.25");
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
    // Run this same suite in both the selected and forced-fallback executables.
    CHECK(std::signbit(parse_number_token<double>("sst", "-0")));
    CHECK(!std::signbit(parse_number_token<double>("sst", "+0")));
    CHECK(parse_number_token<double>("sst", "0e9999") == 0);
    CHECK(parse_number_token<double>("sst", "0e-9999") == 0);
    for (const auto& sample : std::vector<std::pair<const char*, double>>{
           {"1e-300", 1e-300}, {"1e300", 1e300},
           {"4.9406564584124654e-324", std::numeric_limits<double>::denorm_min()},
           {"2.2250738585072014e-308", std::numeric_limits<double>::min()},
           {"1.7976931348623157e308", std::numeric_limits<double>::max()}}) {
      errno = EDOM;
      CHECK(parse_number_token<double>("sst", sample.first) == sample.second);
      CHECK(errno == EDOM);
    }
    for (const char* token : {"0x1p2", "1e+", ".", "nan", "inf", "1e999", "1e-999", "1,5"}) {
      errno = EDOM;
      rejects([&] { parse_number_token<double>("sst", token); });
      CHECK(errno == EDOM);
    }
    rejects([&] { parse_number_token<double>("sst", std::string_view("1\0.5", 4)); });
    const int rounding = std::fegetround();
    for (const int mode : {FE_DOWNWARD, FE_UPWARD}) {
      CHECK(std::fesetround(mode) == 0);
      rejects([&] { parse_number_token<double>("sst", "1e999"); });
      rejects([&] { parse_number_token<double>("sst", "-1e999"); });
    }
    CHECK(std::fesetround(rounding) == 0);
    // Use a non-C process locale when installed, without requiring one on CI.
    const std::string c_locale = std::setlocale(LC_NUMERIC, nullptr);
    bool tested_c_locale = false;
    for (const char* name : {"de_DE.UTF-8", "de_DE.utf8", "German_Germany.1252"}) {
      if (!std::setlocale(LC_NUMERIC, name)) continue;
      CHECK(parse_number_token<double>("sst", "1.5") == 1.5);
      rejects([&] { parse_number_token<double>("sst", "1,5"); });
      tested_c_locale = true;
      break;
    }
    CHECK(std::setlocale(LC_NUMERIC, c_locale.c_str()) != nullptr);
    if (!tested_c_locale) std::cout << "Non-C process locale unavailable; locale-specific check skipped\n";
    std::cout << "AVS numeric strings: mixed separators, strict tokens, bounds and locale passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
