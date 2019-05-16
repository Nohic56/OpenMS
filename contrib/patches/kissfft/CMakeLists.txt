cmake_minimum_required(VERSION 2.8.0)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS ON)

project(kissfft C)

set(KISSFFT_FULL_VERSION "130")

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build a shared library form of kissfft" ON)
endif()

if(NOT WIN32)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -O3 -mtune=native -ffast-math -fomit-frame-pointer -Dkiss_fft_scalar=double -funroll-loops -dA -fverbose-asm")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR} "tools")


# create debug libs with d suffix
set(CMAKE_DEBUG_POSTFIX d)

set(kissfft_lib_sources
  kiss_fft.c
  tools/kfc.c
  tools/kiss_fftnd.c
  tools/kiss_fftndr.c
  tools/kiss_fftr.c
)



set(kissfft_public_headers
  kiss_fft.h
  tools/kfc.h
  tools/kiss_fftnd.h
  tools/kiss_fftndr.h
  tools/kiss_fftr.h
)

set(kissfft_private_headers
  _kiss_fft_guts.h
)

add_library(libkfft ${kissfft_lib_sources} ${kissfft_public_headers} ${kissfft_private_headers})

# embed version in dll
if(NOT CYGWIN)
    set_target_properties(libkfft PROPERTIES VERSION ${KISSFFT_FULL_VERSION})
endif()

# remove lib prefix as the target is already named kissfft
if(NOT MSVC)
   set_target_properties(libkfft PROPERTIES PREFIX "")
endif()

## executables
target_link_libraries(libkfft)

# the install target
install(TARGETS libkfft
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib )

install(FILES ${kissfft_public_headers} DESTINATION include)

# if MSVC we want the dll to be also in the lib-folder
if (MSVC AND BUILD_SHARED_LIBS)
	install(TARGETS libkfft DESTINATION lib	)
endif (MSVC AND BUILD_SHARED_LIBS)
