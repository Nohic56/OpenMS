##################################################
###       GLPK     														 ###
##################################################

MACRO( OPENMS_CONTRIB_BUILD_GLPK )
  OPENMS_LOGHEADER_LIBRARY("glpk")
  #extract: (takes very long.. so skip if possible)
  if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif()
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_GLPK "GLPK" "README")

  ## build the obj/lib
  if (MSVC)
    message(STATUS "Building GLPK library (nmake) ... ")

		configure_file("${PROJECT_SOURCE_DIR}/patches/gplk/Win_config_VC" "${GLPK_DIR}/w32/config.h" COPYONLY)
		configure_file("${PROJECT_SOURCE_DIR}/patches/gplk/Makefile_VC" "${GLPK_DIR}/w32/Makefile_VC" COPYONLY)
		
	  exec_program(${NMAKE_EXECUTABLE} "${GLPK_DIR}/w32"  ## w32 and w64 contain the same makefiles, so it does not matter (the cmd environment we are using to call nmake matters though)
	    ARGS /f Makefile_VC DEBUG=1 check
      OUTPUT_VARIABLE GLPK_MAKE_OUT
	    RETURN_VALUE BUILD_SUCCESS)

    # logfile
    file(APPEND ${LOGFILE} ${GLPK_MAKE_OUT})

	  if (NOT BUILD_SUCCESS EQUAL 0)
      message(STATUS "Building GLPK library (nmake) .. failed")
	    message(FATAL_ERROR ${GLPK_MAKE_OUT})
    else()
      message(STATUS "Building GLPK library (nmake) .. done")
	  endif()

	  exec_program(${NMAKE_EXECUTABLE} "${GLPK_DIR}/w32"  ## w32 and w64 contain the same makefiles, so it does not matter (the cmd environment we are using to call nmake matters though)
	    ARGS /f Makefile_VC check
      OUTPUT_VARIABLE GLPK_MAKE_OUT
	    RETURN_VALUE BUILD_SUCCESS)

    # logfile
    file(APPEND ${LOGFILE} ${GLPK_MAKE_OUT})

	  if (NOT BUILD_SUCCESS EQUAL 0)
      message(STATUS "Building GLPK library (nmake) .. failed")
	    message(FATAL_ERROR ${GLPK_MAKE_OUT})
    else()
      message(STATUS "Building GLPK library (nmake) .. done")
	  endif()

		configure_file("${GLPK_DIR}/src/glpk.h" "${CMAKE_BINARY_DIR}/include/glpk.h" COPYONLY)

  else()  
		# CFLAGS for libsvm compiler (see libsvm Makefile)
	  set(GLPK_CFLAGS "-Wall -O3 -fPIC")
	  
		# add OS X specific flags
		if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
			set(GLPK_CFLAGS "${GLPK_CFLAGS} ${OSX_LIB_FLAG} ${OSX_DEPLOYMENT_FLAG} ${OSX_SYSROOT_FLAG}")
		endif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )

		if (BUILD_SHARED_LIBRARIES)
			set(STATIC_BUILD "--enable-static=no")
			set(SHARED_BUILD "--enable-shared=yes")
		else()
			set(STATIC_BUILD "--enable-static=yes")
			set(SHARED_BUILD "--enable-shared=no")		
		endif()
	

		message( STATUS "Configuring GLPK (./configure --disable-dependency-tracking --prefix ${CMAKE_BINARY_DIR} --disable-shared CPPFLAGS='${GLPK_CFLAGS}' CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} .. ")
    exec_program("./configure" "${GLPK_DIR}"
      ARGS
      --prefix ${CMAKE_BINARY_DIR}
      ${STATIC_BUILD}
      ${SHARED_BUILD}
			--disable-dependency-tracking
			CPPFLAGS='${GLPK_CFLAGS}'
      CXX=${CMAKE_CXX_COMPILER}
      CC=${CMAKE_C_COMPILER}
      OUTPUT_VARIABLE GLPK_CONFIGURE_OUT
      RETURN_VALUE GLPK_CONFIGURE_SUCCESS
      )

    # logfile
    file(APPEND ${LOGFILE} ${GLPK_CONFIGURE_OUT})

    if( NOT GLPK_CONFIGURE_SUCCESS EQUAL 0)
      message( STATUS "Configuring GLPK .. failed")
      message( FATAL_ERROR ${GLPK_CONFIGURE_OUT})
    else()
      message( STATUS "Configuring GLPK .. done")
    endif()

    # make 
    message( STATUS "Building GLPK library (make).. ")
    exec_program(${CMAKE_MAKE_PROGRAM} "${GLPK_DIR}"
      ARGS -j ${NUMBER_OF_JOBS}
      OUTPUT_VARIABLE GLPK_MAKE_OUT
      RETURN_VALUE GLPK_MAKE_SUCCESS
      )

    # logfile
    file(APPEND ${LOGFILE} ${GLPK_MAKE_OUT})

    if( NOT GLPK_MAKE_SUCCESS EQUAL 0)
      message( STATUS "Building GLPK library (make) .. failed")
      message( FATAL_ERROR ${GLPK_MAKE_OUT})
    else()
      message( STATUS "Building GLPK library (make) .. done")
    endif()
    
    # make install
    message( STATUS "Installing GLPK library (make install) .. ")
    exec_program(${CMAKE_MAKE_PROGRAM} "${GLPK_DIR}"
      ARGS "install"
      -j ${NUMBER_OF_JOBS}
      OUTPUT_VARIABLE GLPK_INSTALL_OUT
      RETURN_VALUE GLPK_INSTALL_SUCCESS
      )
    
    # logfile
    file(APPEND ${LOGFILE} ${GLPK_INSTALL_OUT})

    if( NOT GLPK_INSTALL_SUCCESS EQUAL 0)
      message( STATUS "Installing GLPK library (make install) .. failed")
      message( FATAL_ERROR ${GLPK_INSTALL_OUT})
    else()
      message( STATUS "Installing GLPK library (make install) .. done")      
    endif()
endif()


ENDMACRO( OPENMS_CONTRIB_BUILD_GLPK )
