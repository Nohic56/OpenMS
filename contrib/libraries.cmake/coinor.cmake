##################################################
###       COIN-OR															 ###
##################################################
## COIN-OR from http://www.coin-or.org/download/source/CoinMP/CoinMP-1.3.3.tgz
## repacked installed files and created VisualStudio 2008 files

MACRO( OPENMS_CONTRIB_BUILD_COINOR)
  OPENMS_LOGHEADER_LIBRARY("COINOR")
  ## extract: (takes very long.. so skip if possible)
  if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif()
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_COINOR "COINOR" "AUTHORS")

  if (MSVC)
		## changes made to COIN-MP solution files (for all 6 libs):
		## - changed name of output lib in debug mode from $(OutDir)\$(ProjectName).lib to  $(OutDir)\$(ProjectName)d.lib
		## - changed used runtime library to dynamic version (Release & Debug mode)
		## - deleted contents of CoinMP-1.3.3\CoinMP\MSVisualStudio\v8\release (there were precompiled dll's and lib's)

    # set(PATCH_FILE "${PATCH_DIR}/coinor/CoinLpIO.diff")
    # set(PATCHED_FILE "${COINOR_DIR}/CoinUtils/src/CoinLpIO.cpp")
    # OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)
	
    set(MSBUILD_ARGS_SLN "${COINOR_DIR}/CoinMP/MSVisualStudio/v${CONTRIB_MSVC_VERSION}/CoinMP.sln")
    set(MSBUILD_ARGS_TARGET "libCbc")
    OPENMS_BUILDLIB("CoinOR-Cbc (Debug)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Debug" COINOR_DIR)
    OPENMS_BUILDLIB("CoinOR-Cbc (Release)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Release" COINOR_DIR)

    set(MSBUILD_ARGS_TARGET "libCgl")
    OPENMS_BUILDLIB("CoinOR-Cgl (Debug)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Debug" COINOR_DIR)
    OPENMS_BUILDLIB("CoinOR-Cgl (Release)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Release" COINOR_DIR)

    set(MSBUILD_ARGS_TARGET "libClp")
    OPENMS_BUILDLIB("CoinOR-Clp (Debug)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Debug" COINOR_DIR)
    OPENMS_BUILDLIB("CoinOR-Clp (Release)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Release" COINOR_DIR)

    set(MSBUILD_ARGS_TARGET "libCoinUtils")
    OPENMS_BUILDLIB("CoinOR-CoinUtils (Debug)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Debug" COINOR_DIR)
    OPENMS_BUILDLIB("CoinOR-CoinUtils (Release)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Release" COINOR_DIR)

    set(MSBUILD_ARGS_TARGET "libOsi")
    OPENMS_BUILDLIB("CoinOR-Osi (Debug)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Debug" COINOR_DIR)
    OPENMS_BUILDLIB("CoinOR-Osi (Release)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Release" COINOR_DIR)

    set(MSBUILD_ARGS_TARGET "libOsiClp")
    OPENMS_BUILDLIB("CoinOR-OsiClp (Debug)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Debug" COINOR_DIR)
    OPENMS_BUILDLIB("CoinOR-OsiClp (Release)" MSBUILD_ARGS_SLN MSBUILD_ARGS_TARGET "Release" COINOR_DIR)
    
    ###################
    ## copy includes ##
    ###################
		Message(STATUS "Copying include files to ./include/coin ... ")
		file(GLOB_RECURSE INC_FILES "${COINOR_DIR}/*.h" "${COINOR_DIR}/*.hpp")
		
		## create the target directory (coin)
		execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_BINARY_DIR}/include/coin/"
											OUTPUT_VARIABLE MAKE_DIR_OUT
											RESULT_VARIABLE MAKE_DIR_SUCCESS)
		if( NOT MAKE_DIR_SUCCESS EQUAL 0)
			message( STATUS "creating ./include/coin .. failed")
			message( FATAL_ERROR ${MAKE_DIR_OUT})
		endif()
		
		## copying
		foreach (FFF ${INC_FILES})
			execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${FFF}" "${PROJECT_BINARY_DIR}/include/coin/"
											OUTPUT_VARIABLE COPY_INC_OUT
											RESULT_VARIABLE COPY_INC_SUCCESS)
			if( NOT COPY_INC_SUCCESS EQUAL 0)
				message( STATUS "Copying ${FFF} to ./include/coin .. failed")
				message( FATAL_ERROR ${COPY_INC_OUT})
			endif()
		endforeach()		

	else()  ## LINUX & Mac
    #set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/CbcEventHandler.hpp.diff")
    #set(PATCHED_FILE "${COINOR_DIR}/Cbc/src/CbcEventHandler.hpp")
    #OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)

    #set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/CoinTypes.hpp.diff")
    #set(PATCHED_FILE "${COINOR_DIR}/CoinUtils/src/CoinTypes.hpp")
    #OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)

    # set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/Makefile.in.diff")
    # set(PATCHED_FILE "${COINOR_DIR}/Makefile.in")
    # OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)

    # set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/MakefileCoinMP.in.diff")
    # set(PATCHED_FILE "${COINOR_DIR}/CoinMP/Makefile.in")
    # OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)  

    set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/CoinUtils.configure.diff")
    set(PATCHED_FILE "${COINOR_DIR}/CoinUtils/configure")
    OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)

    set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/Clp.configure.diff")
    set(PATCHED_FILE "${COINOR_DIR}/Clp/configure")
    OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)

    set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/Cgl.configure.diff")
    set(PATCHED_FILE "${COINOR_DIR}/Cgl/configure")
    OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)

    set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/CoinMP.configure.diff")
    set(PATCHED_FILE "${COINOR_DIR}/CoinMP/configure")
    OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)

    set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/Cbc.configure.diff")
    set(PATCHED_FILE "${COINOR_DIR}/Cbc/configure")
    OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)

    set(PATCH_FILE "${PROJECT_SOURCE_DIR}/patches/coinor/Osi.configure.diff")
    set(PATCHED_FILE "${COINOR_DIR}/Osi/configure")
    OPENMS_PATCH( PATCH_FILE COINOR_DIR PATCHED_FILE)
  
    # configure -- 
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
      set(COINOR_EXTRA_FLAGS "ADD_FFLAGS='${OSX_DEPLOYMENT_FLAG}' ADD_CFLAGS='${OSX_DEPLOYMENT_FLAG} ${OSX_SYSROOT_FLAG}' ADD_CXXFLAGS='${OSX_LIB_FLAG} ${OSX_DEPLOYMENT_FLAG} ${OSX_SYSROOT_FLAG} -fPIC' --disable-dependency-tracking")
    else()
      set(COINOR_EXTRA_FLAGS "ADD_CXXFLAGS='-fPIC'")
    endif()    

		# check if we prefer shared or static libs
		if (BUILD_SHARED_LIBRARIES)
			set(STATIC_BUILD "--enable-static=no")
			set(SHARED_BUILD "--enable-shared=yes")
		else()
			set(STATIC_BUILD "--enable-static=yes")
			set(SHARED_BUILD "--enable-shared=no")		
		endif()
		
    message( STATUS "Configure COIN-OR library (./configure -C --prefix=${PROJECT_BINARY_DIR} ${STATIC_BUILD} ${SHARED_BUILD} --with-lapack=no --with-blas=no ${COINOR_EXTRA_FLAGS} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER})")
    exec_program("./configure" "${COINOR_DIR}"
      ARGS 
      -C 
      --prefix=${PROJECT_BINARY_DIR}
      ## Following two lines can be combined with prefix
      ## But maybe they avoid building the doc into share (wanted?)
      #--libdir=${CONTRIB_BIN_LIB_DIR} 
      #--includedir=${CONTRIB_BIN_INCLUDE_DIR}
      ${STATIC_BUILD}
      ${SHARED_BUILD}
      --with-lapack=no
      --with-blas=no
      ${COINOR_EXTRA_FLAGS} 
      CXX=${CMAKE_CXX_COMPILER} 
      CC=${CMAKE_C_COMPILER}
      OUTPUT_VARIABLE COINOR_CONFIGURE_OUT
      RETURN_VALUE COINOR_CONFIGURE_SUCCESS
      )
    
    ## logfile
    file(APPEND ${LOGFILE} ${COINOR_CONFIGURE_OUT})

    if( NOT COINOR_CONFIGURE_SUCCESS EQUAL 0)
      message( STATUS "Configure COIN-OR library (./configure -C --prefix=${PROJECT_BINARY_DIR} ${STATIC_BUILD} ${SHARED_BUILD} --with-lapack=no --with-blas=no ${COINOR_EXTRA_FLAGS} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}) .. failed")
      message( FATAL_ERROR ${COINOR_CONFIGURE_OUT})
    else()
      message( STATUS "Configure COIN-OR library (./configure -C --prefix=${PROJECT_BINARY_DIR} ${STATIC_BUILD} ${SHARED_BUILD} --with-lapack=no --with-blas=no ${COINOR_EXTRA_FLAGS} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}) .. done")
    endif()

    ## make install
    message( STATUS "Building and installing COIN-OR library (make install).. ")
    exec_program(${CMAKE_MAKE_PROGRAM} "${COINOR_DIR}"
      ARGS 
      install
      #-j ${NUMBER_OF_JOBS} # the project has problems with multiple jobs. It tries to create folders at the same time and fails.
      -j 1
      OUTPUT_VARIABLE COINOR_MAKE_OUT
      RETURN_VALUE COINOR_MAKE_SUCCESS
      )

    ## logfile
    file(APPEND ${LOGFILE} ${COINOR_MAKE_OUT})

    if( NOT COINOR_MAKE_SUCCESS EQUAL 0)
      message( STATUS "Building and installing COIN-OR library (make install) .. failed")
      message( FATAL_ERROR ${COINOR_MAKE_OUT})
    else()
      message( STATUS "Building and installing COIN-OR library (make install) .. done")
    endif()
  endif()

ENDMACRO( OPENMS_CONTRIB_BUILD_COINOR )
