##################################################
###       Wildmagic   												 ###
##################################################

## build and install eigen
macro( OPENMS_CONTRIB_BUILD_WILDMAGIC )

  OPENMS_LOGHEADER_LIBRARY("wildmagic")

	if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif()
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_WILDMAGIC "WILDMAGIC" "LICENSE_1_0.txt")

  # patch files => remove wildmagic asserts and fix c'tor warning in Vector
  set(_PATCH_FILE "${PATCH_DIR}/wildmagic/Wm5Assert.h.patch")
  set(_PATCHED_FILE "${WILDMAGIC_DIR}/LibCore/Assert/Wm5Assert.h")
  OPENMS_PATCH( _PATCH_FILE WILDMAGIC_DIR _PATCHED_FILE)

  set(_PATCH_FILE "${PATCH_DIR}/wildmagic/Wm5Assert.cpp.patch")
  set(_PATCHED_FILE "${WILDMAGIC_DIR}/LibCore/Assert/Wm5Assert.cpp")
  OPENMS_PATCH( _PATCH_FILE WILDMAGIC_DIR _PATCHED_FILE)

  set(_PATCH_FILE "${PATCH_DIR}/wildmagic/Wm5Vector2.patch")
  set(_PATCHED_FILE "${WILDMAGIC_DIR}/LibMathematics/Algebra/Wm5Vector2.inl")
  OPENMS_PATCH( _PATCH_FILE WILDMAGIC_DIR _PATCHED_FILE)

  set(_PATCH_FILE "${PATCH_DIR}/wildmagic/Wm5CoreLIB.patch")
  set(_PATCHED_FILE "${WILDMAGIC_DIR}/LibCore/Wm5CoreLIB.h")
  OPENMS_PATCH( _PATCH_FILE WILDMAGIC_DIR _PATCHED_FILE)

  # copy build-system addon
	configure_file(${PROJECT_SOURCE_DIR}/patches/wildmagic/root_CMakeLists.txt ${WILDMAGIC_DIR}/CMakeLists.txt COPYONLY)
	configure_file(${PROJECT_SOURCE_DIR}/patches/wildmagic/libcore_CMakeLists.txt ${WILDMAGIC_DIR}/LibCore/CMakeLists.txt COPYONLY)
	configure_file(${PROJECT_SOURCE_DIR}/patches/wildmagic/libmath_CMakeLists.txt ${WILDMAGIC_DIR}/LibMathematics/CMakeLists.txt COPYONLY)


  # we want out-of-source builds
  set(_WILDMAGIC_BUILD_DIR "${WILDMAGIC_DIR}/build")
  file(TO_NATIVE_PATH "${_WILDMAGIC_BUILD_DIR}" _WILDMAGIC_NATIVE_BUILD_DIR)

  if(APPLE)
    set( _WILDMAGIC_CMAKE_ARGS
    	"-D CMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}"
    	"-D CMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}"
    	"-D CMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}"
    )
  else()
    set( _WILDMAGIC_CMAKE_ARGS "")
  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${_WILDMAGIC_NATIVE_BUILD_DIR})

	message(STATUS "Generating wildmagic build system .. ")
  execute_process(COMMAND ${CMAKE_COMMAND}
												-G "${CMAKE_GENERATOR}"
                        -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
												-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                        ${_WILDMAGIC_CMAKE_ARGS}
												${WILDMAGIC_DIR}
									WORKING_DIRECTORY ${_WILDMAGIC_NATIVE_BUILD_DIR}
									OUTPUT_VARIABLE _WILDMAGIC_CMAKE_OUT
									ERROR_VARIABLE _WILDMAGIC_CMAKE_ERR
									RESULT_VARIABLE _WILDMAGIC_CMAKE_SUCCESS)

  # output to logfile
  file(APPEND ${LOGFILE} ${_WILDMAGIC_CMAKE_OUT})
  file(APPEND ${LOGFILE} ${_WILDMAGIC_CMAKE_ERR})

	if (NOT _WILDMAGIC_CMAKE_SUCCESS EQUAL 0)
		message(FATAL_ERROR "Generating wildmagic build system .. failed")
	else()
		message(STATUS "Generating wildmagic build system .. done")
	endif()

  # the install target on windows has a different name then on mac/lnx
  if(MSVC)
      set(_WILDMAGIC_INSTALL_TARGET "INSTALL")
  else()
      set(_WILDMAGIC_INSTALL_TARGET "install")
  endif()

  # build release first
	message(STATUS "Building wildmagic library (Release) .. ")
	execute_process(COMMAND ${CMAKE_COMMAND} --build ${_WILDMAGIC_NATIVE_BUILD_DIR} --target ${_WILDMAGIC_INSTALL_TARGET} --config Release
									WORKING_DIRECTORY ${_WILDMAGIC_NATIVE_BUILD_DIR}
									OUTPUT_VARIABLE _WILDMAGIC_BUILD_OUT
									ERROR_VARIABLE _WILDMAGIC_BUILD_ERR
									RESULT_VARIABLE _WILDMAGIC_BUILD_SUCCESS)

	# output to logfile
	file(APPEND ${LOGFILE} ${_WILDMAGIC_BUILD_OUT})
	file(APPEND ${LOGFILE} ${_WILDMAGIC_BUILD_ERR})

	if (NOT _WILDMAGIC_BUILD_SUCCESS EQUAL 0)
		message(FATAL_ERROR "Building wildmagic library (Release) .. failed")
	else()
		message(STATUS "Building wildmagic library (Release) .. done")
	endif()

  # we also want the debug lib on windows
  if(MSVC)
    # build debug
  	message(STATUS "Building wildmagic library (Debug) .. ")
  	execute_process(COMMAND ${CMAKE_COMMAND} --build ${_WILDMAGIC_NATIVE_BUILD_DIR} --target ${_WILDMAGIC_INSTALL_TARGET} --config Debug
  									WORKING_DIRECTORY ${_WILDMAGIC_NATIVE_BUILD_DIR}
  									OUTPUT_VARIABLE _WILDMAGIC_BUILD_OUT
  									ERROR_VARIABLE _WILDMAGIC_BUILD_ERR
  									RESULT_VARIABLE _WILDMAGIC_BUILD_SUCCESS)

  	# output to logfile
  	file(APPEND ${LOGFILE} ${_WILDMAGIC_BUILD_OUT})
  	file(APPEND ${LOGFILE} ${_WILDMAGIC_BUILD_ERR})

  	if (NOT _WILDMAGIC_BUILD_SUCCESS EQUAL 0)
  		message(FATAL_ERROR "Building wildmagic library (Debug) .. failed")
  	else()
  		message(STATUS "Building wildmagic library (Debug) .. done")
  	endif()
  endif()
endmacro( OPENMS_CONTRIB_BUILD_WILDMAGIC )