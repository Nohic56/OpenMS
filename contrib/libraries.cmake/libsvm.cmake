##################################################
###       LIBSVM   														 ###
##################################################

MACRO( OPENMS_CONTRIB_BUILD_LIBSVM )
  OPENMS_LOGHEADER_LIBRARY("libSVM")
  #extract: (takes very long.. so skip if possible)
  if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif()
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_LIBSVM "LIBSVM" "README")

	## we use our own CMakeLists.txt for libsvn
	## if bzip is update, ensure that the CMakeLists.txt is still valid
	configure_file(${PROJECT_SOURCE_DIR}/patches/libsvm/CMakeLists.txt ${LIBSVM_DIR}/CMakeLists.txt COPYONLY)

  #patch...
  if (MSVC)
	  set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/libsvm/svm.cpp.diff_vs")
	  set(PATCHED_FILE "${LIBSVM_DIR}/svm.cpp")
	  OPENMS_PATCH( PATCH_FILE LIBSVM_DIR PATCHED_FILE)
  elseif (MINGW)
	  set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/libsvm/svm.cpp.diff_mingw")
	  set(PATCHED_FILE "${LIBSVM_DIR}/svm.cpp")
	  OPENMS_PATCH( PATCH_FILE LIBSVM_DIR PATCHED_FILE)
  endif (MSVC)

  ## build the obj/lib
  if (MSVC)
		message(STATUS "Generating libsvm build system .. ")
    execute_process(COMMAND ${CMAKE_COMMAND}
													-D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
													-G "${CMAKE_GENERATOR}"
													-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
													.
										WORKING_DIRECTORY ${LIBSVM_DIR}
										OUTPUT_VARIABLE LIBSVM_CMAKE_OUT
										ERROR_VARIABLE LIBSVM_CMAKE_ERR
										RESULT_VARIABLE LIBSVM_CMAKE_SUCCESS)

		# output to logfile
		file(APPEND ${LOGFILE} ${LIBSVM_CMAKE_OUT})
		file(APPEND ${LOGFILE} ${LIBSVM_CMAKE_ERR})

		if (NOT LIBSVM_CMAKE_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Generating libsvm build system .. failed")
		else()
			message(STATUS "Generating libsvm build system .. done")
		endif()

		message(STATUS "Building libsvm lib (Debug) .. ")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${LIBSVM_DIR} --target INSTALL --config Debug
										WORKING_DIRECTORY ${LIBSVM_DIR}
										OUTPUT_VARIABLE LIBSVM_BUILD_OUT
										ERROR_VARIABLE LIBSVM_BUILD_ERR
										RESULT_VARIABLE LIBSVM_BUILD_SUCCESS)

		# output to logfile
		file(APPEND ${LOGFILE} ${LIBSVM_BUILD_OUT})
		file(APPEND ${LOGFILE} ${LIBSVM_BUILD_ERR})

		if (NOT LIBSVM_BUILD_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Building libsvm lib (Debug) .. failed")
		else()
			message(STATUS "Building libsvm lib (Debug) .. done")
		endif()

		## rebuild as release
		message(STATUS "Building libsvm lib (Release) .. ")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${LIBSVM_DIR} --target INSTALL --config Release
										WORKING_DIRECTORY ${LIBSVM_DIR}
										OUTPUT_VARIABLE LIBSVM_BUILD_OUT
										ERROR_VARIABLE LIBSVM_BUILD_ERR
										RESULT_VARIABLE LIBSVM_BUILD_SUCCESS)
		# output to logfile
		file(APPEND ${LOGFILE} ${LIBSVM_BUILD_OUT})
		file(APPEND ${LOGFILE} ${LIBSVM_BUILD_ERR})

		if (NOT LIBSVM_BUILD_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Building libsvm lib (Release) .. failed")
		else()
			message(STATUS "Building libsvm lib (Release) .. done")
		endif()
  else()  
		message(STATUS "Generating libsvm build system .. ")
		if (APPLE)
 			execute_process(COMMAND ${CMAKE_COMMAND}
														-D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
														-G "${CMAKE_GENERATOR}"
														-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
														-D CMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
														-D CMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
														-D CMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
														.
											WORKING_DIRECTORY ${LIBSVM_DIR}
											OUTPUT_VARIABLE LIBSVM_CMAKE_OUT
											ERROR_VARIABLE LIBSVM_CMAKE_ERR
											RESULT_VARIABLE LIBSVM_CMAKE_SUCCESS)
		else()
			execute_process(COMMAND ${CMAKE_COMMAND}
														-D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
														-G "${CMAKE_GENERATOR}"
														-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
														.
											WORKING_DIRECTORY ${LIBSVM_DIR}
											OUTPUT_VARIABLE LIBSVM_CMAKE_OUT
											ERROR_VARIABLE LIBSVM_CMAKE_ERR
											RESULT_VARIABLE LIBSVM_CMAKE_SUCCESS)
		endif()
		
		# output to logfile
		file(APPEND ${LOGFILE} ${LIBSVM_CMAKE_OUT})
		file(APPEND ${LOGFILE} ${LIBSVM_CMAKE_ERR})

    if (NOT LIBSVM_CMAKE_SUCCESS EQUAL 0)
      message(FATAL_ERROR "Generating libsvm build system .. failed")
    else()
      message(STATUS "Generating libsvm build system .. done")
    endif()
    
		## rebuild as release
		message(STATUS "Building libsvm lib (Release) .. ")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${LIBSVM_DIR} --target install --config Release
										WORKING_DIRECTORY ${LIBSVM_DIR}
										OUTPUT_VARIABLE LIBSVM_BUILD_OUT
										ERROR_VARIABLE LIBSVM_BUILD_ERR
										RESULT_VARIABLE LIBSVM_BUILD_SUCCESS)
		# output to logfile
		file(APPEND ${LOGFILE} ${LIBSVM_BUILD_OUT})
		file(APPEND ${LOGFILE} ${LIBSVM_BUILD_ERR})

		if (NOT LIBSVM_BUILD_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Building libsvm lib (Release) .. failed")
		else()
			message(STATUS "Building libsvm lib (Release) .. done")
		endif()
endif()

## copy includes
configure_file(${LIBSVM_DIR}/svm.h ${CONTRIB_BIN_INCLUDE_DIR} COPYONLY)

ENDMACRO( OPENMS_CONTRIB_BUILD_LIBSVM )
