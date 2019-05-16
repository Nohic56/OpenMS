##################################################
###       BZIP2   														 ###
##################################################

MACRO( OPENMS_CONTRIB_BUILD_BZIP2 )
  OPENMS_LOGHEADER_LIBRARY("bzip2")

	if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif()
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_BZIP2 "bzip2" "README")

	## we use our own CMakeLists.txt for bzip2
	## if bzip is update, ensure that the CMakeLists.txt is still valid
	configure_file(${PROJECT_SOURCE_DIR}/patches/bzip2/CMakeLists.txt ${BZIP2_DIR}/CMakeLists.txt COPYONLY)

  ## build the obj/lib
  if (MSVC)
		message(STATUS "Generating bzip2 build system .. ")
    execute_process(COMMAND ${CMAKE_COMMAND}
													-D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
													-G "${CMAKE_GENERATOR}"
													-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
													.
										WORKING_DIRECTORY ${BZIP2_DIR}
										OUTPUT_VARIABLE BZIP_CMAKE_OUT
										ERROR_VARIABLE BZIP_CMAKE_ERR
										RESULT_VARIABLE BZIP_CMAKE_SUCCESS)

		# output to logfile
		file(APPEND ${LOGFILE} ${BZIP_CMAKE_OUT})
		file(APPEND ${LOGFILE} ${BZIP_CMAKE_ERR})

		if (NOT BZIP_CMAKE_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Generating bzip2 build system .. failed")
		else()
			message(STATUS "Generating bzip2 build system .. done")
		endif()

		message(STATUS "Building bzip2 lib (Debug) .. ")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${BZIP2_DIR} --target INSTALL --config Debug
										WORKING_DIRECTORY ${BZIP2_DIR}
										OUTPUT_VARIABLE BZIP_BUILD_OUT
										ERROR_VARIABLE BZIP_BUILD_ERR
										RESULT_VARIABLE BZIP_BUILD_SUCCESS)

		# output to logfile
		file(APPEND ${LOGFILE} ${BZIP_BUILD_OUT})
		file(APPEND ${LOGFILE} ${BZIP_BUILD_ERR})

		if (NOT BZIP_BUILD_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Building bzip2 lib (Debug) .. failed")
		else()
			message(STATUS "Building bzip2 lib (Debug) .. done")
		endif()

		## rebuild as release
		message(STATUS "Building bzip2 lib (Release) .. ")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${BZIP2_DIR} --target INSTALL --config Release
										WORKING_DIRECTORY ${BZIP2_DIR}
										OUTPUT_VARIABLE BZIP_BUILD_OUT
										ERROR_VARIABLE BZIP_BUILD_ERR
										RESULT_VARIABLE BZIP_BUILD_SUCCESS)
		# output to logfile
		file(APPEND ${LOGFILE} ${BZIP_BUILD_OUT})
		file(APPEND ${LOGFILE} ${BZIP_BUILD_ERR})

		if (NOT BZIP_BUILD_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Building bzip2 lib (Release) .. failed")
		else()
			message(STATUS "Building bzip2 lib (Release) .. done")
		endif()

  else()

    # CFLAGS for libsvm compiler (see libsvm Makefile)
    set(BZIP2_CFLAGS "-Wall -Winline -O2 -g -fPIC -D_FILE_OFFSET_BITS=64")

    # add OS X specific flags
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
      set(BZIP2_CFLAGS "${BZIP2_CFLAGS} ${OSX_DEPLOYMENT_FLAG}")
    endif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )

		message(STATUS "Generating bzip2 build system .. ")
		if (APPLE)
			execute_process(COMMAND ${CMAKE_COMMAND}
														-D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
														-D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
														-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
														-D CMAKE_C_FLAGS='${BZIP2_CFLAGS}'
														-G "${CMAKE_GENERATOR}"
														-D CMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
														-D CMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
														-D CMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
														.
											WORKING_DIRECTORY ${BZIP2_DIR}
											OUTPUT_VARIABLE BZIP_CMAKE_OUT
											ERROR_VARIABLE BZIP_CMAKE_ERR
											RESULT_VARIABLE BZIP_CMAKE_SUCCESS)
		else()
			execute_process(COMMAND ${CMAKE_COMMAND}
												-D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
												-D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
												-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
												-D CMAKE_C_FLAGS='${BZIP2_CFLAGS}'
												-G "${CMAKE_GENERATOR}"
												.
									WORKING_DIRECTORY ${BZIP2_DIR}
									OUTPUT_VARIABLE BZIP_CMAKE_OUT
									ERROR_VARIABLE BZIP_CMAKE_ERR
									RESULT_VARIABLE BZIP_CMAKE_SUCCESS)
		endif()
		# output to logfile
		file(APPEND ${LOGFILE} ${BZIP_CMAKE_OUT})
		file(APPEND ${LOGFILE} ${BZIP_CMAKE_ERR})

		if (NOT BZIP_CMAKE_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Generating bzip2 build system .. failed")
		else()
			message(STATUS "Generating bzip2 build system .. done")
		endif()

		message(STATUS "Building bzip2 lib (Release) .. ")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${BZIP2_DIR} --target install --config Release
										WORKING_DIRECTORY ${BZIP2_DIR}
										OUTPUT_VARIABLE BZIP_BUILD_OUT
										ERROR_VARIABLE BZIP_BUILD_ERR
										RESULT_VARIABLE BZIP_BUILD_SUCCESS)
		# output to logfile
		file(APPEND ${LOGFILE} ${BZIP_BUILD_OUT})
		file(APPEND ${LOGFILE} ${BZIP_BUILD_ERR})

		if (NOT BZIP_BUILD_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Building bzip2 lib (Release) .. failed")
		else()
			message(STATUS "Building bzip2 lib (Release) .. done")
		endif()
endif()

ENDMACRO( OPENMS_CONTRIB_BUILD_BZIP2 )
