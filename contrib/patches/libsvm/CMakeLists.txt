cmake_minimum_required(VERSION 2.8.0)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS ON)

project(libsvm)

set(LIBSVM_FULL_VERSION "3.12")

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build a shared library." ON)
endif()

if(WIN32)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS} -nologo -MD -O2 -EHsc -I. -D __WIN32__ -D _CRT_SECURE_NO_DEPRECATE")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS} -nologo -MDd -Od -EHsc -I. -D __WIN32__ -D _CRT_SECURE_NO_DEPRECATE")	
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -O3 -fPIC")
endif()

# create debug libs with d suffix
set(CMAKE_DEBUG_POSTFIX d)

# the libsvm sources
set(libsvm_sources
	svm.cpp
)

# add def file for dll
if (WIN32 AND BUILD_SHARED_LIBS)
	set(libsvm_sources 
    ${libsvm_sources} 
    svm.def
  )
endif()

# the libsvm header
set(libsvm_public_headers
  svm.h
)

## libsvm library
add_library(libsvm ${libsvm_sources})

# set version for embedding in the SO/DYLIB
set_target_properties(libsvm PROPERTIES 
	VERSION "${LIBSVM_FULL_VERSION}"
	SOVERSION "${LIBSVM_FULL_VERSION}")

# remove default prefix of cmake, since our target is already named libsvm
set_target_properties(libsvm PROPERTIES PREFIX "")

# add the libsvm apps
set(svm_executables
	svm-predict
	svm-train
	svm-scale
)

foreach(i ${svm_executables})
  add_executable(${i} ${i}.c)
  target_link_libraries(${i} libsvm)
endforeach(i)

# the install target
install(TARGETS libsvm ${svm_executables}
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib )

install(FILES ${libsvm_public_headers} DESTINATION include)
	
# if MSVC we want the dll to be also in the lib-folder
# for linking against libsvm
if (MSVC AND BUILD_SHARED_LIBS)
	install(TARGETS libsvm DESTINATION lib	)
endif (MSVC AND BUILD_SHARED_LIBS)
