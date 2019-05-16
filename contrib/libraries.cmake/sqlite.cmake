##################################################
###       SQLITE                               ###
##################################################

MACRO( OPENMS_CONTRIB_BUILD_SQLITE )
  OPENMS_LOGHEADER_LIBRARY("SQLITE")
  #extract: (takes very long.. so skip if possible)
  if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif(MSVC)
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_SQLITE "SQLITE" "INSTALL")
  
  if(MSVC)
    message( STATUS "Building SQLITE library in  ${SQLITE_DIR}")
	
	if (${BUILD_SHARED_LIBRARIES})
		# there is a Makefile.am, but it uses (broken) .def export instead of just a simple call to cl. So we just use:
		execute_process(COMMAND "cl" "sqlite3.c" "-DSQLITE_API=__declspec(dllexport)" "-link" "-dll" "-out:sqlite3.dll"
		               WORKING_DIRECTORY "${SQLITE_DIR}"
		               RESULT_VARIABLE _SQLITE_RES
		               OUTPUT_VARIABLE _SQLITE_OUT
		               ERROR_VARIABLE _SQLITE_ERR
		               )
					   
	else()
		execute_process(COMMAND "cl" "/c" "sqlite3.c"
						WORKING_DIRECTORY "${SQLITE_DIR}"
						RESULT_VARIABLE _SQLITE_RES
						OUTPUT_VARIABLE _SQLITE_OUT
						ERROR_VARIABLE _SQLITE_ERR
						)
		execute_process(COMMAND "lib" "sqlite3.obj"
						WORKING_DIRECTORY "${SQLITE_DIR}"
						RESULT_VARIABLE _SQLITE_RES
						OUTPUT_VARIABLE _SQLITE_OUT
						ERROR_VARIABLE _SQLITE_ERR
						)
	endif()
	
    if (NOT _SQLITE_RES EQUAL 0)
      message( STATUS "Building sqlite failed")
      file(APPEND ${LOGFILE} "sqlite failed" )
	  message( FATAL_ERROR ${_SQLITE_OUT})
    else()
      message( STATUS "Building sqlite worked")
      file(APPEND ${LOGFILE} "sqlite worked" )
    endif()
    
    file(APPEND ${LOGFILE} ${_SQLITE_ERR})
    file(APPEND ${LOGFILE} ${_SQLITE_OUT})
    
    configure_file(${SQLITE_DIR}/sqlite3.h ${PROJECT_BINARY_DIR}/include/sqlite/sqlite3.h COPYONLY)
	if(${BUILD_SHARED_LIBRARIES})
		configure_file(${SQLITE_DIR}/sqlite3.dll ${PROJECT_BINARY_DIR}/lib/sqlite.dll COPYONLY)
	endif()
    configure_file(${SQLITE_DIR}/sqlite3.lib ${PROJECT_BINARY_DIR}/lib/sqlite.lib COPYONLY)
	
  else()
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
     set (SQLITE_CUSTOM_FLAGS "${CXX_OSX_FLAGS}")
    endif()
    
    # configure -- 
    set( ENV{CC} ${CMAKE_C_COMPILER} )
    set( ENV{CXX} ${CMAKE_CXX_COMPILER} )
    set( ENV{CFLAGS} ${SQLITE_CUSTOM_FLAGS})

    if (BUILD_SHARED_LIBRARIES)
      set(STATIC_BUILD "--enable-static=no")
      set(SHARED_BUILD "--enable-shared=yes")
    else()
      set(STATIC_BUILD "--enable-static=yes")
      set(SHARED_BUILD "--enable-shared=no")        
    endif()
        
    message( STATUS "Configure SQLITE library (./configure --prefix ${CMAKE_BINARY_DIR} --with-pic --disable-shared) .. ")
    exec_program("sed" "${SQLITE_DIR}"
        ARGS
        -i "\"s/PACKAGE_STRING='sqlite 3.15.0'/PACKAGE_STRING='sqlite3.15.0'/\"" configure
        OUTPUT_VARIABLE SQLITE_SED_OUT
        RETURN_VALUE SQLITE_SED_SUCCESS
    )

    exec_program("./configure" "${SQLITE_DIR}"
      ARGS
      --prefix ${CMAKE_BINARY_DIR}
      --with-pic
      ${STATIC_BUILD}
      ${SHARED_BUILD}
      OUTPUT_VARIABLE SQLITE_CONFIGURE_OUT
      RETURN_VALUE SQLITE_CONFIGURE_SUCCESS
      )

    # logfile
    file(APPEND ${LOGFILE} ${SQLITE_CONFIGURE_OUT})

    if( NOT SQLITE_CONFIGURE_SUCCESS EQUAL 0)
      message( STATUS "Configure SQLITE library (./configure --prefix ${CMAKE_BINARY_DIR} --with-pic --disable-shared) .. failed")
      message( FATAL_ERROR ${SQLITE_CONFIGURE_OUT})
    else()
      message( STATUS "Configure SQLITE library (./configure --prefix ${CMAKE_BINARY_DIR} --with-pic --disable-shared) .. done")
    endif()
  
    # make 
    message( STATUS "Building SQLITE library (make) .. ")
    exec_program(${CMAKE_MAKE_PROGRAM} "${SQLITE_DIR}"
      ARGS -j ${NUMBER_OF_JOBS}
      OUTPUT_VARIABLE SQLITE_MAKE_OUT
      RETURN_VALUE SQLITE_MAKE_SUCCESS
      )

    file(APPEND ${LOGFILE} ${SQLITE_MAKE_OUT})

    if( NOT SQLITE_MAKE_SUCCESS EQUAL 0)
      message( STATUS "Building SQLITE library (make) .. failed")
      message( FATAL_ERROR ${SQLITE_MAKE_OUT})
    else()
      message( STATUS "Building SQLITE library (make) .. done")
    endif()

    # make install
    message( STATUS "Installing SQLITE library (make install) .. ")
    exec_program(${CMAKE_MAKE_PROGRAM} "${SQLITE_DIR}"
      ARGS "install"
      -j ${NUMBER_OF_JOBS}
      OUTPUT_VARIABLE SQLITE_INSTALL_OUT
      RETURN_VALUE SQLITE_INSTALL_SUCCESS
      )

    file(APPEND ${LOGFILE} ${SQLITE_INSTALL_OUT})

    if( NOT SQLITE_INSTALL_SUCCESS EQUAL 0)
      message( STATUS "Installing SQLITE library (make install) .. failed")      
      message( FATAL_ERROR ${SQLITE_INSTALL_OUT})
    else()
      message( STATUS "Installing SQLITE library (make install) .. done")
    endif()
endif()

ENDMACRO( OPENMS_CONTRIB_BUILD_SQLITE )
