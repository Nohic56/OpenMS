##################################################
###       XERCES 															 ###
##################################################
## XERCES from http://xerces.apache.org/xerces-c/

MACRO( OPENMS_CONTRIB_BUILD_XERCESC )
	OPENMS_LOGHEADER_LIBRARY("Xerces-C")
	#extract: (takes very long.. so skip if possible)
	if(MSVC)
		set(ZIP_ARGS "x -y -osrc")
	else()
		set(ZIP_ARGS "xzf")
	endif()
	OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_XERCES "XERCES" "CREDITS")
	## Regarding transcoder choices. Let Xerces figure it out.

	if (WIN32)
		#set(XERCES_EXTRA_CMAKE_FLAGS "-D...")
		message( STATUS "Generating XERCES-C cmake build system for Debug and Release...")
		execute_process(COMMAND ${CMAKE_COMMAND}
								-D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
								-D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
								-D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
								-G "${CMAKE_GENERATOR}"
								-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
								-D CMAKE_POSITION_INDEPENDENT_CODE=ON
								#${XERCES_EXTRA_CMAKE_FLAGS}
								.
						WORKING_DIRECTORY ${XERCES_DIR}
						OUTPUT_VARIABLE XERCES_CMAKE_OUT
						ERROR_VARIABLE XERCES_CMAKE_ERR
						RESULT_VARIABLE XERCES_CMAKE_SUCCESS)				
						
		# logfile
		file(APPEND ${LOGFILE} ${XERCES_CMAKE_OUT} ${XERCES_Release_ERR})
		
		execute_process(COMMAND ${CMAKE_COMMAND}
						--build .
						--config Debug
						--target install
				WORKING_DIRECTORY ${XERCES_DIR}
				OUTPUT_VARIABLE XERCES_Debug_OUT
				ERROR_VARIABLE XERCES_Debug_ERR
				RESULT_VARIABLE XERCES_Debug_SUCCESS)
				
		file(APPEND ${LOGFILE} ${XERCES_Debug_OUT})
				
		execute_process(COMMAND ${CMAKE_COMMAND}
						--build .
						--config Release
						--target install
				WORKING_DIRECTORY ${XERCES_DIR}
				OUTPUT_VARIABLE XERCES_Release_OUT
				ERROR_VARIABLE XERCES_Release_ERR
				RESULT_VARIABLE XERCES_Release_SUCCESS)
		
		file(APPEND ${LOGFILE} ${XERCES_Release_OUT} ${XERCES_Release_ERR})
	else()
		## Linux and Mac
                set(_PATCH_FILE "${PATCH_DIR}/xercesc/xerces_dynlib_link_suffix.patch")
		set(_PATCHED_FILE "${XERCES_DIR}/src/CMakeLists.txt")
		OPENMS_PATCH( _PATCH_FILE XERCES_DIR _PATCHED_FILE)
	        
		set(_transcoder "gnuiconv") # to use on linux to reduce external libs (even if ICU is found)
		
		if(APPLE)
	            set(XERCES_EXTRA_CMAKE_FLAGS "-DCMAKE_SHARED_LIBARY_SUFFIX=dylib")
		    set(_transcoder "macosunicodeconverter") # standard on macOS
                endif()
		
		## Release
		message( STATUS "Reconfiguring XERCES-C cmake build system for Release...")
		execute_process(COMMAND ${CMAKE_COMMAND}
                                                                -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                                                -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
								-D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBRARIES}
								-G "${CMAKE_GENERATOR}"
								-D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
								-D CMAKE_INSTALL_LIBDIR=lib
								-D CMAKE_BUILD_TYPE=Release
								-D CMAKE_POSITION_INDEPENDENT_CODE=ON
								-D CMAKE_CXX_FLAGS="-fPIC"
								-Dnetwork:BOOL=OFF
								-Dtranscoder=${_transcoder}
								${XERCES_EXTRA_CMAKE_FLAGS}
								.
						WORKING_DIRECTORY ${XERCES_DIR}
						OUTPUT_VARIABLE XERCES_CMAKE_Release_OUT
						ERROR_VARIABLE XERCES_CMAKE_Release_ERR
						RESULT_VARIABLE XERCES_CMAKE_Release_SUCCESS)
						
		file(APPEND ${LOGFILE} ${XERCES_CMAKE_Release_OUT} ${XERCES_Release_ERR})
						
		message( STATUS "Building and installing XERCES-C for Release...")				
		execute_process(COMMAND ${CMAKE_COMMAND}
								--build .
								--target install
						WORKING_DIRECTORY ${XERCES_DIR}
						OUTPUT_VARIABLE XERCES_Release_OUT
						ERROR_VARIABLE XERCES_Release_ERR
						RESULT_VARIABLE XERCES_Release_SUCCESS)				

		file(APPEND ${LOGFILE} ${XERCES_Release_OUT} ${XERCES_Release_ERR})

	endif()
ENDMACRO( OPENMS_CONTRIB_BUILD_XERCESC )
