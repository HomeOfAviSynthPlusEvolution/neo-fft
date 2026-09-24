include(CheckCXXSourceCompiles)

option(NEO_FFT_FORCE_DOUBLE_PARSE_FALLBACK "Use C-locale strtod instead of double from_chars" OFF)

# Check declarations and symbols against the selected SDK/deployment target.
# A toolchain may otherwise turn try_compile into a compile-only static library.
function(neo_fft_check_double_parsing)
  set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)
  check_cxx_source_compiles("
    #include <charconv>
    int main(int argc, char** argv) {
      double value = 0;
      const char* text = argv[argc - 1];
      auto result = std::from_chars(text, text + 1, value, std::chars_format::general);
      return int(result.ec) + (value == 1.5);
    }
  " NEO_FFT_HAS_DOUBLE_FROM_CHARS)
  if(BUILD_TESTING OR NEO_FFT_FORCE_DOUBLE_PARSE_FALLBACK OR NOT NEO_FFT_HAS_DOUBLE_FROM_CHARS)
    if(WIN32)
      set(locale_probe "
        #include <cstdlib>
        #include <locale.h>
        int main(int argc, char** argv) {
          auto locale = _create_locale(LC_NUMERIC, \"C\");
          char* end;
          double value = _strtod_l(argv[argc - 1], &end, locale);
          _free_locale(locale);
          return value == 1.5;
        }")
    else()
      set(locale_probe "
        #include <cstdlib>
        #include <locale.h>
        #if defined(__APPLE__) || defined(__FreeBSD__)
        #include <xlocale.h>
        #endif
        int main(int argc, char** argv) {
          auto locale = newlocale(LC_NUMERIC_MASK, \"C\", nullptr);
          char* end;
          double value = strtod_l(argv[argc - 1], &end, locale);
          freelocale(locale);
          return value == 1.5;
        }")
    endif()
    check_cxx_source_compiles("${locale_probe}" NEO_FFT_HAS_LOCALE_STRTOD)
    if(NOT NEO_FFT_HAS_LOCALE_STRTOD AND
       (NEO_FFT_FORCE_DOUBLE_PARSE_FALLBACK OR NOT NEO_FFT_HAS_DOUBLE_FROM_CHARS))
      message(FATAL_ERROR "AviSynth numeric strings need double std::from_chars or C-locale strtod_l")
    endif()
  endif()
endfunction()

neo_fft_check_double_parsing()

function(neo_fft_configure_double_parsing target)
  if(NEO_FFT_HAS_DOUBLE_FROM_CHARS AND NOT NEO_FFT_FORCE_DOUBLE_PARSE_FALLBACK)
    target_compile_definitions(${target} PRIVATE NEO_FFT_HAS_DOUBLE_FROM_CHARS=1)
  endif()
endfunction()
