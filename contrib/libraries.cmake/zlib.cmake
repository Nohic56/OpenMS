##################################################
###       ZLIB   							   ###
##################################################

MACRO( OPENMS_CONTRIB_BUILD_ZLIB )
  OPENMS_LOGHEADER_LIBRARY("zlib")
  #extract: (takes very long.. so skip if possible)
  if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif()
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_ZLIB "ZLIB" "README")
	
  ## build the obj/lib
  if (MSVC)		
    message(STATUS "Generating zlib build system .. ")
    execute_process(COMMAND ${CMAKE_COMMAND}
                            -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
                            -D INSTALL_BIN_DIR=${PROJECT_BINARY_DIR}/lib
                            -G "${CMAKE_GENERATOR}"
                            -D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                            ${ZLIB_EXTRA_CMAKE_FLAG}
                            .
                    WORKING_DIRECTORY ${ZLIB_DIR}
                    OUTPUT_VARIABLE ZLIB_CMAKE_OUT
                    ERROR_VARIABLE ZLIB_CMAKE_ERR
                    RESULT_VARIABLE ZLIB_CMAKE_SUCCESS)

		# output to logfile
		file(APPEND ${LOGFILE} ${ZLIB_CMAKE_OUT})
		file(APPEND ${LOGFILE} ${ZLIB_CMAKE_ERR})

		if (NOT ZLIB_CMAKE_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Generating zlib build system .. failed")
		else()
			message(STATUS "Generating zlib build system .. done")
		endif()

		message(STATUS "Building zlib lib (Debug) .. ")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${ZLIB_DIR} --target INSTALL --config Debug
										WORKING_DIRECTORY ${ZLIB_DIR}
										OUTPUT_VARIABLE ZLIB_BUILD_OUT
										ERROR_VARIABLE ZLIB_BUILD_ERR
										RESULT_VARIABLE ZLIB_BUILD_SUCCESS)

		# output to logfile
		file(APPEND ${LOGFILE} ${ZLIB_BUILD_OUT})
		file(APPEND ${LOGFILE} ${ZLIB_BUILD_ERR})

		if (NOT ZLIB_BUILD_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Building zlib lib (Debug) .. failed")
		else()
			message(STATUS "Building zlib lib (Debug) .. done")
		endif()

		# rebuild as release
		message(STATUS "Building zlib lib (Release) .. ")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${ZLIB_DIR} --target INSTALL --config Release
										WORKING_DIRECTORY ${ZLIB_DIR}
										OUTPUT_VARIABLE ZLIB_BUILD_OUT
										ERROR_VARIABLE ZLIB_BUILD_ERR
										RESULT_VARIABLE ZLIB_BUILD_SUCCESS)
		# output to logfile
		file(APPEND ${LOGFILE} ${ZLIB_BUILD_OUT})
		file(APPEND ${LOGFILE} ${ZLIB_BUILD_ERR})

		if (NOT ZLIB_BUILD_SUCCESS EQUAL 0)
			message(FATAL_ERROR "Building zlib lib (Release) .. failed")
		else()
			message(STATUS "Building zlib lib (Release) .. done")
		endif()
  
	else() ## Linux/MacOS  

    # CFLAGS for libsvm compiler (see libsvm Makefile)
    set(ZLIB_CFLAGS "-Wall" "-O3" "-fPIC")

    message(STATUS "Generating zlib build system .. ")
    execute_process(COMMAND ${CMAKE_COMMAND}
                            -G "${CMAKE_GENERATOR}"
                            -D CMAKE_BUILD_TYPE=Release
                            -D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                            -D CMAKE_C_FLAGS=${ZLIB_CFLAGS}
                            .
                    WORKING_DIRECTORY ${ZLIB_DIR}
                    OUTPUT_VARIABLE ZLIB_CMAKE_OUT
                    ERROR_VARIABLE ZLIB_CMAKE_ERR
                    RESULT_VARIABLE ZLIB_CMAKE_SUCCESS)

    # rebuild as release
    message(STATUS "Building zlib lib (Release) .. ")
    execute_process(COMMAND ${CMAKE_COMMAND} --build ${ZLIB_DIR} --target install
                    WORKING_DIRECTORY ${ZLIB_DIR}
                    OUTPUT_VARIABLE ZLIB_BUILD_OUT
                    ERROR_VARIABLE ZLIB_BUILD_ERR
                    RESULT_VARIABLE ZLIB_BUILD_SUCCESS)
    # output to logfile
    file(APPEND ${LOGFILE} ${ZLIB_BUILD_OUT})
    file(APPEND ${LOGFILE} ${ZLIB_BUILD_ERR})

    if(APPLE) ## Somehow libz does not add the path to the install_name. -> Problems while resolving dependecies during shipping
    execute_process(COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id ${PROJECT_BINARY_DIR}/lib/libz.1.dylib ${PROJECT_BINARY_DIR}/lib/libz.1.dylib)
    endif()

    if (NOT ZLIB_BUILD_SUCCESS EQUAL 0)
      message(FATAL_ERROR "Building zlib lib (Release) .. failed")
    else()
      message(STATUS "Building zlib lib (Release) .. done")
    endif()

endif()

ENDMACRO( OPENMS_CONTRIB_BUILD_ZLIB )
