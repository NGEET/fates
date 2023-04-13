string(APPEND CFLAGS " -qopt-report -march=core-avx2 -mtune=core-avx2 -no-fma")
string(APPEND FFLAGS " -qopt-report -march=core-avx2 -mtune=core-avx2 -no-fma")
if (DEBUG)
  string(APPEND CMAKE_OPTS " -DPIO_ENABLE_LOGGING=ON")
endif()
set(SLIBS "-lnetcdf -lnetcdff -lpnetcdf -lblas -llapack -lpthread -lm -ldl")
string(APPEND CPPDEFS " -DNO_SHR_VMATH")

