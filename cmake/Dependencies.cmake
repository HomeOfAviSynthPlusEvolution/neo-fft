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
  GIT_TAG f1d51bd0217d3878f995375e95c2827b4604facf
  PATCH_COMMAND "${CMAKE_COMMAND}" "-DNEO_FFT_DUALSYNTH_SOURCE=<SOURCE_DIR>"
    -P "${CMAKE_CURRENT_LIST_DIR}/PatchDualSynth.cmake"
  SOURCE_SUBDIR "${_ds_subdir}")
FetchContent_Declare(pocketfft
  GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
  GIT_TAG c90e55b3d529f8efa40ed01a20de22405f45fc65
  SOURCE_SUBDIR neo_fft_header_only)
FetchContent_MakeAvailable(dualsynth2 pocketfft)
unset(_ds_subdir)
