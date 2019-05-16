##################################################
###       BOOST                                ###
##################################################

MACRO( OPENMS_CONTRIB_BUILD_BOOST)
  OPENMS_LOGHEADER_LIBRARY("BOOST")
  
  set( BOOST_BUILD_TYPE "static")
  if (BUILD_SHARED_LIBRARIES)
    set( BOOST_BUILD_TYPE "shared")
  endif()
  
  ## extract boost library
  if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif()
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_BOOST "BOOST" "index.htm")
  
  if(MSVC) ## build boost library for windows
    
    ## omitting the version (i.e. 'toolset=msvc'), causes Boost to use the latest(!) VS it can find the system -- irrespective of the current env (and its cl.exe)
    set(TOOLSET "toolset=msvc-${CONTRIB_MSVC_VERSION}.0") 
    
    if (NOT QUICKBUILD)
      ## not a Visual Studio project .. just build by hand
      message(STATUS "Bootstrapping Boost libraries (bootstrap.bat) ...")
      execute_process(COMMAND bootstrap.bat
                      WORKING_DIRECTORY ${BOOST_DIR}
                      OUTPUT_VARIABLE BOOST_BOOTSTRAP_OUT
                      ERROR_VARIABLE BOOST_BOOTSTRAP_OUT   # use same variable for stderr as stdout to merge streams
                      RESULT_VARIABLE BOOST_BOOTSTRAP_SUCCESS)
      
      file(APPEND  ${LOGFILE} ${BOOST_BOOTSTRAP_OUT})
      
      ## check for failed bootstrapping. Even if failing the return code can be 0 (indicating success), so we additionally check the output 
      if ((NOT BOOST_BOOTSTRAP_SUCCESS EQUAL 0) OR (BOOST_BOOTSTRAP_OUT MATCHES "[fF]ailed"))
        message(STATUS "Bootstrapping Boost libraries (bootstrap.bat) ... failed\nOutput was:\n ${BOOST_BOOTSTRAP_OUT}\nEnd of output.\n")
        message(STATUS "Renaming bootstrap.log to bootstrap_firstTry.log")
        file(RENAME ${BOOST_DIR}/bootstrap.log ${BOOST_DIR}/bootstrap_firstTry.log)
        ### on some command lines bootstrapping fail (e.g. the toolset is too new) or will give:
        # "Building Boost.Build engine. The input line is too long."
        ## ,thus we provide a backup bjam.exe(32bit), which hopefully works on all target systems.
        ## However this bjam results in a version mismatch and a warning (that you can ignore).
        message(STATUS " ... trying fallback with backup bjam.exe ...")
        configure_file("${PROJECT_SOURCE_DIR}/patches/boost/bjam.exe" "${BOOST_DIR}/bjam.exe" COPYONLY)
      else()
        message(STATUS "Bootstrapping Boost libraries (bootstrap.bat) ... done")
      endif()

      set(BOOST_CMD_ARGS "${BOOST_ARG}" 
                         "install" 
                         "-j${NUMBER_OF_JOBS}" 
                         "--prefix=${PROJECT_BINARY_DIR}" 
                         "--layout=tagged"                   # create libnames without vcXXX in filename; include dir is /include/boost (as opposed to "versioned" where /include/boost-1.52/boost plus ...vc110.lib
                         "--with-math" 
                         "--with-date_time" 
                         "--with-iostreams" 
                         "--with-regex"
                         "--with-system"
                         "--with-thread"
                         "--build-type=complete"
                         "--disable-icu"
                         "-sZLIB_SOURCE=${ZLIB_DIR}"
                         "-sBZIP2_SOURCE=${BZIP2_DIR}" 
                         "runtime-link=shared"
                         "link=${BOOST_BUILD_TYPE}" 
                         "${TOOLSET}")

     ## WARNING: boost call is not "space in path" save yet (the easy way of using \" does not work out of the box
     message(STATUS "Building Boost library (bjam ${BOOST_CMD_ARGS}) .. ")
     execute_process(COMMAND bjam.exe ${BOOST_CMD_ARGS}
                     WORKING_DIRECTORY ${BOOST_DIR}
                     OUTPUT_VARIABLE BUILD_BOOST_OUT
                     ERROR_VARIABLE BUILD_BOOST_ERR
                     RESULT_VARIABLE BUILD_BOOST)
      
     # output to logfile
     file(APPEND ${LOGFILE} ${BUILD_BOOST_OUT})

     if (NOT BUILD_BOOST EQUAL 0)
       message(STATUS "Building Boost library (bjam ${BOOST_CMD_ARGS}) .. failed")
       message(STATUS "BUILD_BOOST_OUT = ${BUILD_BOOST_OUT}")
       message(STATUS "BUILD_BOOST_ERR = ${BUILD_BOOST_ERR}")
       message(STATUS "BUILD_BOOST = ${BUILD_BOOST}")
       message(FATAL_ERROR ${BUILD_BOOST_OUT})
     else()
       message(STATUS "Building Boost library (bjam ${BOOST_CMD_ARGS}) .. done")
     endif()
    
   endif() ## end quickbuild
 
  else() ## LINUX/MAC

    # we need to know the compiler version for proper formating boost user-config.jam
    determine_compiler_version()

    # use proper toolchain (random guesses. There is not proper documentation)
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
      if(APPLE)
        set(_boost_bootstrap_toolchain "darwin")
        set(_boost_toolchain "clang-darwin")
      else()
        set(_boost_bootstrap_toolchain "clang")
        set(_boost_toolchain "clang")
      endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      if(APPLE)
        ## For Apples old GCC (tag in lib name will be xgcc)
        set(_boost_bootstrap_toolchain "darwin")
        set(_boost_toolchain "darwin") 
      else()
        set(_boost_bootstrap_toolchain "gcc")
        set(_boost_toolchain "gcc")
      endif()
    endif()

    ## In case we specified a non-existent toolchain. Just use whatever CMake detected/knows.
    ##file(APPEND ${BOOST_DIR}/tools/build/src/user-config.jam
    ##  "using ${_boost_toolchain} : ${CXX_COMPILER_VERSION_MAJOR}.${CXX_COMPILER_VERSION_MINOR} : ${CMAKE_CXX_COMPILER};\n")

    if(APPLE AND CMAKE_OSX_DEPLOYMENT_TARGET)
      ## Boost looks for installed SDKs, but sometimes you dont have them. Add them still to not fail. Clang will handle it.
      file(APPEND ${BOOST_DIR}/tools/build/src/tools/darwin.jam
        "feature.extend macosx-version-min : ${CMAKE_OSX_DEPLOYMENT_TARGET} ;\n")

      ## Add corresponding linker flags. Empty is not possible, therefore the if.
      if(OSX_LIB_FLAG)
        set(BOOST_LINKER_FLAGS linkflags=${OSX_LIB_FLAG})
      endif()
    endif()
    

    # bootstrap boost
    message(STATUS "Bootstrapping Boost libraries (./bootstrap.sh --prefix=${PROJECT_BINARY_DIR} --with-toolset=${_boost_bootstrap_toolchain} --with-libraries=date_time,iostreams,math,regex,system,thread) ...")
    execute_process(COMMAND ./bootstrap.sh --prefix=${PROJECT_BINARY_DIR} --with-toolset=${_boost_bootstrap_toolchain} --with-libraries=iostreams,math,date_time,regex,system,thread
                    WORKING_DIRECTORY ${BOOST_DIR}
                    OUTPUT_VARIABLE BOOST_BOOTSTRAP_OUT
                    ERROR_VARIABLE BOOST_BOOTSTRAP_OUT
                    RESULT_VARIABLE BOOST_BOOTSTRAP_SUCCESS)

    # logfile
    file(APPEND ${LOGFILE} ${BOOST_BOOTSTRAP_OUT})
    if (NOT BOOST_BOOTSTRAP_SUCCESS EQUAL 0)
      message(STATUS "Bootstrapping Boost libraries (./bootstrap.sh --prefix=${PROJECT_BINARY_DIR} --with-libraries=iostreams,math,date_time,regex,system,thread) ... failed")
      message(FATAL_ERROR ${BOOST_BOOTSTRAPPING_OUT})
    else()
      message(STATUS "Bootstrapping Boost libraries (./bootstrap.sh --prefix=${PROJECT_BINARY_DIR} --with-libraries=iostreams,math,date_time,regex,system,thread) ... done")
    endif()

    # boost cmd
    set (BOOST_CMD "./bjam toolset=${_boost_toolchain} -j ${NUMBER_OF_JOBS} --disable-icu -sZLIB_SOURCE=${ZLIB_DIR} -sBZIP2_SOURCE=${BZIP2_DIR} ${OSX_DEPLOYMENT_TARGET_STRING} link=${BOOST_BUILD_TYPE} cxxflags='-fPIC ${OSX_LIB_FLAGS}' ${BOOST_LINKER_FLAGS} install --build-type=complete --layout=tagged --threading=single,multi")
    
    # boost install
    message(STATUS "Installing Boost libraries (${BOOST_CMD}) ...")
    execute_process(COMMAND ./bjam toolset=${_boost_toolchain} -j ${NUMBER_OF_JOBS} --disable-icu -sZLIB_SOURCE="${ZLIB_DIR}" -sBZIP2_SOURCE="${BZIP2_DIR}" ${OSX_DEPLOYMENT_TARGET_STRING} link=${BOOST_BUILD_TYPE} cxxflags=-fPIC ${OSX_CXX_FLAG} ${BOOST_LINKER_FLAGS} install --build-type=complete --layout=tagged --threading=single,multi
                    WORKING_DIRECTORY ${BOOST_DIR}
                    OUTPUT_VARIABLE BOOST_INSTALL_OUT
                    ERROR_VARIABLE BOOST_INSTALL_OUT
                    RESULT_VARIABLE BOOST_INSTALL_SUCCESS)

    # logfile 
    file(APPEND ${LOGFILE} ${BOOST_INSTALL_OUT})
    if (NOT BOOST_INSTALL_SUCCESS EQUAL 0)
      message(STATUS "Installing Boost libraries (${BOOST_CMD}) ... failed")
      message(FATAL_ERROR ${BOOST_INSTALL_OUT})
    else()
      message(STATUS "Installing Boost libraries (${BOOST_CMD}) ... done")
    endif()
  endif()

ENDMACRO(OPENMS_CONTRIB_BUILD_BOOST)
