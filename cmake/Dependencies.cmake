include(FetchContent)
if(NEO_FFT_BUILD_VAPOURSYNTH OR NEO_FFT_BUILD_AVISYNTH)
  set(DS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(DS_BUILD_ACCEPTANCE_PLUGIN OFF CACHE BOOL "" FORCE)
  set(DS_ENABLE_AVISYNTH ${NEO_FFT_BUILD_AVISYNTH} CACHE BOOL "" FORCE)
  set(DS_ENABLE_VAPOURSYNTH ${NEO_FFT_BUILD_VAPOURSYNTH} CACHE BOOL "" FORCE)
  set(_ds_subdir .)
else()
  set(_ds_subdir neo_fft_header_only)
endif()
FetchContent_Declare(dualsynth2
  GIT_REPOSITORY https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2.git
  GIT_TAG 627127332a7a360f8cbf644be545115bc16585cd
  SOURCE_SUBDIR "${_ds_subdir}")
FetchContent_Declare(pocketfft
  GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
  GIT_TAG c90e55b3d529f8efa40ed01a20de22405f45fc65
  SOURCE_SUBDIR neo_fft_header_only)
FetchContent_MakeAvailable(dualsynth2 pocketfft)
unset(_ds_subdir)
