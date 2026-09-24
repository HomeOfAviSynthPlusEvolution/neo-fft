# The pinned DualSynth2 revision converts negative row offsets to size_t.
# Patch only FetchContent's checkout; caller-supplied source overrides are owned
# by the caller and must provide their own equivalent negative-stride fix.
set(_header "${NEO_FFT_DUALSYNTH_SOURCE}/include/dualsynth/span2d.hpp")
file(READ "${_header}" _original)
set(_patched "${_original}")

function(replace_checked old new)
  string(FIND "${_patched}" "${old}" _old_pos)
  string(FIND "${_patched}" "${new}" _new_pos)
  if(NOT _old_pos EQUAL -1 AND _new_pos EQUAL -1)
    string(REPLACE "${old}" "${new}" _patched "${_patched}")
    set(_patched "${_patched}" PARENT_SCOPE)
  elseif(_old_pos EQUAL -1 AND NOT _new_pos EQUAL -1)
    # Already patched: FetchContent can repeat the patch step.
  else()
    message(FATAL_ERROR "Unexpected DualSynth2 span2d.hpp; recheck the negative-stride patch")
  endif()
endfunction()

replace_checked(
  "return data_[static_cast<std::size_t>(y) * stride_ + x];"
  "return data_[static_cast<std::ptrdiff_t>(y) * stride_ + static_cast<std::ptrdiff_t>(x)];")
replace_checked(
  "return data_ + static_cast<std::size_t>(y) * stride_;"
  "return data_ + static_cast<std::ptrdiff_t>(y) * stride_;")
replace_checked(
  "return BasicPlane(data_ + static_cast<std::size_t>(y) * stride_ + x, w, h, stride_bytes());"
  "return BasicPlane(data_ + static_cast<std::ptrdiff_t>(y) * stride_ + x, w, h, stride_bytes());")

if(NOT _patched STREQUAL _original)
  file(WRITE "${_header}" "${_patched}")
endif()
