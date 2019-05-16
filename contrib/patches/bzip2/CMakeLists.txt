cmake_minimum_required(VERSION 2.8.0)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS ON)

project(bzip2 C)

set(BZIP2_FULL_VERSION "1.0.5")

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build a shared library form of zlib" ON)
endif()

if(NOT WIN32)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Winline -O2 -g -D_FILE_OFFSET_BITS=64 ")
endif()
	
if(MSVC)
	add_definitions(/D_FILE_OFFSET_BITS=64x)
endif(MSVC)
	
# create debug libs with d suffix
set(CMAKE_DEBUG_POSTFIX d)

set(bzip2_lib_sources 
  # 	@cat words0
  blocksort.c
  huffman.c
  crctable.c
  randtable.c
  compress.c
  decompress.c
  bzlib.c
)

# we need the def file on windows for correct exports
if(MSVC AND BUILD_SHARED_LIBS)
  set(bzip2_lib_sources 
    ${bzip2_lib_sources} 
    libbz2.def
  )
endif(MSVC AND BUILD_SHARED_LIBS)

set(bzip2_public_headers
  bzlib.h
)

set(bzip2_private_headers
  bzlib_private.h
)

add_library(libbz2 ${bzip2_lib_sources} ${bzip2_public_headers} ${bzip2_private_headers})

# embed version in dll 
if(NOT CYGWIN)
    set_target_properties(libbz2 PROPERTIES VERSION ${BZIP2_FULL_VERSION})
endif()

# remove lib prefix as the target is already named libbz2
if(NOT MSVC)
   set_target_properties(libbz2 PROPERTIES PREFIX "")
endif()

## executables
add_executable(bzip2 bzip2.c)
target_link_libraries(bzip2 libbz2)

add_executable(bzip2recover bzip2recover.c)

# the install target
install(TARGETS bzip2 bzip2recover libbz2
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib )

install(FILES ${bzip2_public_headers} DESTINATION include)
	
# if MSVC we want the dll to be also in the lib-folder
# for linking against libsvm
if (MSVC AND BUILD_SHARED_LIBS)
	install(TARGETS libbz2 DESTINATION lib	)
endif (MSVC AND BUILD_SHARED_LIBS)
