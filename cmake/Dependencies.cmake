include(FetchContent)
if(NEO_FFT_BUILD_VAPOURSYNTH)
  set(DS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(DS_BUILD_ACCEPTANCE_PLUGIN OFF CACHE BOOL "" FORCE)
  set(DS_ENABLE_AVISYNTH OFF CACHE BOOL "" FORCE)
  set(DS_ENABLE_VAPOURSYNTH ON CACHE BOOL "" FORCE)
  set(_ds_subdir .)
else()
  set(_ds_subdir neo_fft_header_only)
endif()
FetchContent_Declare(dualsynth2
  GIT_REPOSITORY https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2.git
  GIT_TAG 2a24d6b4fe808692bfa10c1f9734a3c50c15774e
  SOURCE_SUBDIR "${_ds_subdir}")
FetchContent_Declare(pocketfft
  GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
  GIT_TAG 5f27d5a8f51c5c25030cb22abf434decc9faf0ff
  SOURCE_SUBDIR neo_fft_header_only)
FetchContent_MakeAvailable(dualsynth2 pocketfft)
unset(_ds_subdir)
