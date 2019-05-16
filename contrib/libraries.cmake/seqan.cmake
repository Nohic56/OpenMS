##################################################
###       SEQAN                                ###
##################################################
##
##

MACRO( OPENMS_CONTRIB_BUILD_SEQAN )
  OPENMS_LOGHEADER_LIBRARY("SeqAn")
  #extract: (takes very long.. so skip if possible)
  if(MSVC)
    set(ZIP_ARGS "x -y -osrc")
  else()
    set(ZIP_ARGS "xzf")
  endif()
  OPENMS_SMARTEXTRACT(ZIP_ARGS ARCHIVE_SEQAN "SEQAN" "include/seqan/version.h")

  # create seqan build system and install to target directory
  set(_SEQAN_BUILD_DIR "${SEQAN_DIR}/build")
  file(TO_NATIVE_PATH "${_SEQAN_BUILD_DIR}" _SEQAN_NATIVE_BUILD_DIR)

  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${_SEQAN_NATIVE_BUILD_DIR})

  message(STATUS "Generating seqan build system .. ")
  execute_process(COMMAND ${CMAKE_COMMAND}
                        -G "${CMAKE_GENERATOR}"
                        -D CMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                        -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY
                        ${SEQAN_DIR}
                  WORKING_DIRECTORY ${_SEQAN_NATIVE_BUILD_DIR}
                  OUTPUT_VARIABLE _SEQAN_CMAKE_OUT
                  ERROR_VARIABLE _SEQAN_CMAKE_ERR
                  RESULT_VARIABLE _SEQAN_CMAKE_SUCCESS)

  # output to logfile
  file(APPEND ${LOGFILE} ${_SEQAN_CMAKE_OUT})
  file(APPEND ${LOGFILE} ${_SEQAN_CMAKE_ERR})

  if (NOT _SEQAN_CMAKE_SUCCESS EQUAL 0)
    message(FATAL_ERROR "Generating seqan build system .. failed")
  else()
    message(STATUS "Generating seqan build system .. done")
  endif()

  # the install target on windows has a different name then on mac/lnx
  if(MSVC)
      set(_SEQAN_INSTALL_TARGET "INSTALL")
  else()
      set(_SEQAN_INSTALL_TARGET "install")
  endif()

  # fake the doc target by simply creating the doc/html folder
  set(_SEQAN_EXPECTED_DOC_DIR "${_SEQAN_BUILD_DIR}/docs/html")
  file(TO_NATIVE_PATH "${_SEQAN_EXPECTED_DOC_DIR}" _SEQAN_NATIVE_EXPECTED_DOC_DIR)
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${_SEQAN_NATIVE_EXPECTED_DOC_DIR})

  message(STATUS "Installing seqan headers .. ")
  execute_process(COMMAND ${CMAKE_COMMAND} --build ${_SEQAN_NATIVE_BUILD_DIR} --target ${_SEQAN_INSTALL_TARGET} --config Release
                  WORKING_DIRECTORY ${_SEQAN_NATIVE_BUILD_DIR}
                  OUTPUT_VARIABLE _SEQAN_BUILD_OUT
                  ERROR_VARIABLE _SEQAN_BUILD_ERR
                  RESULT_VARIABLE _SEQAN_BUILD_SUCCESS)

  # output to logfile
  file(APPEND ${LOGFILE} ${_SEQAN_BUILD_OUT})
  file(APPEND ${LOGFILE} ${_SEQAN_BUILD_ERR})

  if (NOT _SEQAN_BUILD_SUCCESS EQUAL 0)
    message(FATAL_ERROR "Installing seqan headers .. failed")
  else()
    message(STATUS "Installing seqan headers .. done")
  endif()

ENDMACRO( OPENMS_CONTRIB_BUILD_SEQAN )
