# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Stephan Aiche, Chris Bielow $
# --------------------------------------------------------------------------

macro(determine_compiler_version )
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR  "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                    OUTPUT_VARIABLE CXX_COMPILER_VERSION)
    string(REGEX MATCHALL "[0-9]+" CXX_COMPILER_COMPONENTS ${CXX_COMPILER_VERSION})
    list(GET CXX_COMPILER_COMPONENTS 0 CXX_COMPILER_VERSION_MAJOR)
    list(LENGTH CXX_COMPILER_COMPONENTS VERSION_COMPONENTS_COUNT)
    # set minor version number of compiler if present, 0 otherwise
    if (VERSION_COMPONENTS_COUNT GREATER 1)
      list(GET CXX_COMPILER_COMPONENTS 1 CXX_COMPILER_VERSION_MINOR)
    else()
      set(CXX_COMPILER_VERSION_MINOR 0)
    endif()
  endif()
endmacro()

## validates the archive for the given library
## @param libname The library that should be validated
macro(validate_archive libname)
  set(_archive_folder "${PROJECT_BINARY_DIR}/archives")
  set(_target_file "${_archive_folder}/${ARCHIVE_${libname}}")
  set(_target_sha1 ${ARCHIVE_${libname}_SHA1})

  message(STATUS "Validating archive for ${libname} .. ")

  file(SHA1 ${_target_file} _downloaded_sha1 )
  if(NOT "${_downloaded_sha1}" STREQUAL "${_target_sha1}")
    file(REMOVE ${_target_file})
    if("${_downloaded_sha1}" STREQUAL "da39a3ee5e6b4b0d3255bfef95601890afd80709")
      message(STATUS "Validating archive for ${libname} .. Found empty archive.")
      message(STATUS "Check your internet connection.")
    else()
      message(STATUS "Validating archive for ${libname} .. sha1 mismatch (expected: ${_target_sha1} got: ${_downloaded_sha1})")
      message(STATUS "The archive file for ${libname} seems to be damaged and will be removed.")
    endif()
    message(STATUS "Please try to rebuild ${libname} to trigger a new download of the archive.")
    message(STATUS "If this fails again, please contact the OpenMS support.")
  else()
    message(STATUS "Validating archive for ${libname} .. successful")
  endif()
endmacro()


## downloads the archive for the given library
## @param libname The library that should be downloaded
macro(download_contrib_archive libname)
  # constant
  # Currently this points to an FTP at FU Berlin
  # Sources are checked out regularly from OpenMS/contrib-sources via a cron job
  set(_BASE_URL "http://ftp.mi.fu-berlin.de/pub/OpenMS/contrib-sources/")

  # the files/folders where downloads are stored
  set(_archive_folder "${PROJECT_BINARY_DIR}/archives")
  set(_target_file "${_archive_folder}/${ARCHIVE_${libname}}")
  set(_target_sha1 ${ARCHIVE_${libname}_SHA1})
  set(_remote_file ${ARCHIVE_${libname}})
  set(_full_url "${_BASE_URL}${_remote_file}")

  # create separate folder for our archives
  if(NOT EXISTS ${_archive_folder})
    file(MAKE_DIRECTORY ${_archive_folder})
  endif(NOT EXISTS ${_archive_folder})

  if(NOT DEFINED RETRIES)
    set(_retries 3)
  elseif(NOT "${RETRIES}" MATCHES "^[1-9][0-9]*")
    message(FATAL_ERROR "RETRIES argument should be a number greater or equal to 1. Found \"${RETRIES}\"")
  else()
    set(_retries ${RETRIES})
  endif()
  
  set(_attempt 0)
  set(_succeeded 0)
  if(NOT EXISTS ${_target_file})
    while(${_attempt} LESS ${_retries} AND NOT ${_succeeded})
      math(EXPR _attempt "${_attempt}+1")
      message(STATUS "Downloading ${libname} .. Try #${_attempt}")
      # download
      file(DOWNLOAD ${_full_url} "${_target_file}")
      message(STATUS "Downloading ${libname} .. done")
      message(STATUS "Validating archive for ${libname} .. ")

      file(SHA1 ${_target_file} _downloaded_sha1 )
      if(NOT "${_downloaded_sha1}" STREQUAL "${_target_sha1}")
        file(REMOVE ${_target_file})
        if("${_downloaded_sha1}" STREQUAL "da39a3ee5e6b4b0d3255bfef95601890afd80709")
          message(STATUS "Validating archive for ${libname} .. Found empty archive.")
        else()
          message(STATUS "Validating archive for ${libname} .. sha1 mismatch (expected: ${_target_sha1} got: ${_downloaded_sha1})")
          message(STATUS "The archive file for ${libname} seems to be damaged and will be removed.")
        endif()
      else()
        message(STATUS "Validating archive for ${libname} .. successful")
        set(_succeeded 1)
      endif()
    endwhile()
    if (NOT ${_succeeded})
      message(FATAL_ERROR "Download and validation failed after ${_attempt} tries. Check your internet connection and try to rebuild ${libname} to trigger a new download of the archive.")
      message(STATUS "If this fails again, please contact the OpenMS support.")
    endif()
  else()
    message(STATUS "Downloading ${libname} .. skipped (already downloaded)")
  endif(NOT EXISTS ${_target_file})
endmacro()


## extract archive to PROJECT_BINARY_DIR
## Warning: it is important that ${${libname}_DIR} exists!
## @param zip_args_varname Name of the variable containing the 7z arguments (without the actual archive)
## @param libfile_varname Name of the variable containing the actual archive
## @param libname Name of the library (should be upper-case -> ${${libname}_DIR} should exist (see top of this file))
## @param checkfile Name of the (last) file to be extracted from the archive (we will use ${${libname}_DIR}/${checkfile} to see if it exists)
MACRO (OPENMS_SMARTEXTRACT zip_args_varname libfile_varname libname checkfile)
  string(TOUPPER ${libname} libnameUP)

  # download the archive
  download_contrib_archive(${libnameUP})

  if(MSVC)
    # ms way
    message(STATUS "Extracting ${libname} ..")
    if (NOT EXISTS ${${libnameUP}_DIR}/${checkfile}) ## last file to be extracted
      exec_program(${PROGRAM_ZIP} ${PROJECT_BINARY_DIR}
        OUTPUT_VARIABLE ZIP_OUT
        ARGS ${${zip_args_varname}} "\"${PROJECT_BINARY_DIR}/archives/${${libfile_varname}}\""
        RETURN_VALUE EXTRACT_SUCCESS)

      if (NOT EXTRACT_SUCCESS EQUAL 0)
        message(STATUS "Extracting ${libname} .. failed")
        message(FATAL_ERROR ${ZIP_OUT})
      elseif(WIN32 AND NOT ${${libfile_varname}_TAR} STREQUAL "") ## maybe we need to extract a tar.gz in two steps
        message(STATUS "Extracting ${libname} .. done (1st pass)")
        message(STATUS "Extracting ${libname} .. ")
        ## extract the tar
        exec_program(${PROGRAM_ZIP} ${PROJECT_BINARY_DIR}
          OUTPUT_VARIABLE ZIP2_OUT
          ARGS ${${zip_args_varname}} "\"${PROJECT_BINARY_DIR}/src/${${libfile_varname}_TAR}\""
          RETURN_VALUE EXTRACT_SUCCESS)

        # logfile
        file(APPEND ${LOGFILE} ${ZIP2_OUT})

        if (NOT EXTRACT_SUCCESS EQUAL 0)
          message(STATUS "Extracting ${libname} .. failed")
          message(FATAL_ERROR ${ZIP2_OUT})
        else()
          message(STATUS "Extracting ${libname} .. done (2nd pass)")
        endif()
      else() ## extraction finished
        message(STATUS "Extracting ${libname} .. done")
      endif()
    else()
      message(STATUS "Extracting ${libname} .. skipped (already exists)")
    endif()
  else()
    ## do it the unix way
    message(STATUS "Extracting ${libname} .. ")
    if (NOT EXISTS ${${libnameUP}_DIR}/${checkfile}) ## last file to be extracted
      exec_program(${PROGRAM_ZIP} ${PROJECT_BINARY_DIR}
        ARGS ${${zip_args_varname}} "${PROJECT_BINARY_DIR}/archives/${${libfile_varname}}" " -C " ${CONTRIB_BIN_SOURCE_DIR}
        OUTPUT_VARIABLE ZIP_OUT
        RETURN_VALUE EXTRACT_SUCCESS)

      if (NOT EXTRACT_SUCCESS EQUAL 0)
        message(STATUS "Extracting ${libname} .. failed")
        message(FATAL_ERROR ${ZIP_OUT})
      else()
        message(STATUS "Extracting ${libname} .. done")
      endif()
    else()
      message(STATUS "Extracting ${libname} .. skipped (already exists)")
    endif()
  endif()
ENDMACRO (OPENMS_SMARTEXTRACT)


## build a library with its VisualStudio solution
## @param libname Name of the library (for text messages only)
## @param solutionfile_varname Name of the variable which hold the path to the solution file (.sln)
## @param target_varname Name of the variable that holds the target to build
## @param config Usually either 'Debug' or 'Release' as string
## @param workingdir_varname Name of the variable that holds the directory where to execute MSBuild in
MACRO ( OPENMS_BUILDLIB libname solutionfile_varname target_varname config workingdir_varname)
  message(STATUS "Building ${libname} ... ")
  file(TO_NATIVE_PATH ${${solutionfile_varname}} _sln_file_path)
  set (MSBUILD_ARGS "/p:Configuration=${config} /consoleloggerparameters:Verbosity=minimal /target:${${target_varname}} /p:Platform=${WIN_PLATFORM_ARG} \"${_sln_file_path}\"")
  if(NOT CONTRIB_MSVC_VERSION STREQUAL "8")
    set (MSBUILD_ARGS "${MSBUILD_ARGS} /maxcpucount")
  endif()

  exec_program(${MSBUILD_EXECUTABLE} ${${workingdir_varname}}
    ARGS ${MSBUILD_ARGS}
    OUTPUT_VARIABLE BUILDLIB_OUT
    RETURN_VALUE BUILD_SUCCESS)

  # logfile
  file(APPEND ${LOGFILE} ${BUILDLIB_OUT})

  if (NOT BUILD_SUCCESS EQUAL 0)
    message(STATUS "Building ${libname} ... failed")
    message(FATAL_ERROR ${BUILDLIB_OUT})
  endif()
ENDMACRO (OPENMS_BUILDLIB)

## patch a file
## @param args_varname Name of the variable that holds all arguments to patch(.exe) [use -i instead of '<' as last argument]
## @param patchfile_varname Name of the variable that holds the patchfile
## @param workingdir_varname Name of the variable that points to the directory where the command is executed
## @param patchedfile_varname Name of the variable that is patched
MACRO ( OPENMS_PATCH patchfile_varname workingdir_varname patchedfile_varname)
  ## First try: with --binary (because of EOL problems, OS and patch.exe specific)
  set( PATCH_ARGUMENTS "-p0 --binary -b -N -i") ## NOTE: always keep -i as last argument !!
  if (EXISTS ${${patchedfile_varname}}.orig)
    message(STATUS "Patching ${${patchedfile_varname}} ... skipped (already applied)")
  else()
    message(STATUS "Try patching ${${patchedfile_varname}} with binary option ... ")
    exec_program(${PROGRAM_PATCH} ${${workingdir_varname}}
      ARGS ${PATCH_ARGUMENTS} "\"${${patchfile_varname}}\""
      OUTPUT_VARIABLE PATCH_OUT
      RETURN_VALUE PATCH_SUCCESS)

    # logfile
    file(APPEND ${LOGFILE} "${PATCH_OUT}\n\r")

    if (NOT PATCH_SUCCESS EQUAL 0)
      ## Second try: without --binary
      set( PATCH_ARGUMENTS "-p0 -b -N -i") ## NOTE: always keep -i as last argument !!
      message(STATUS "Try patching ${${patchedfile_varname}} without binary option ... ")
      exec_program(${PROGRAM_PATCH} ${${workingdir_varname}}
        ARGS ${PATCH_ARGUMENTS} "\"${${patchfile_varname}}\""
        OUTPUT_VARIABLE PATCH_OUT
        RETURN_VALUE PATCH_SUCCESS)

      # logfile
      file(APPEND ${LOGFILE} "${PATCH_OUT}\n\r")
    endif()
        
    if (NOT PATCH_SUCCESS EQUAL 0)
      message(STATUS "Patching ${${patchedfile_varname}} ... failed (with and without --binary option)")
      message(STATUS "Check if the patch was created with 'diff -u' and if the paths are correct!")
      message(STATUS "Call was: ${${workingdir_varname}}: ${PROGRAM_PATCH}  ${PATCH_ARGUMENTS} \"${${patchfile_varname}}\"")
      message(FATAL_ERROR ${PATCH_OUT})
    else()
      message(STATUS "Patching ${${patchedfile_varname}} ... done")
    endif()
  endif()
ENDMACRO (OPENMS_PATCH)


## copy a directory
## @param dir_source Name of the variable that holds the source directory
## @param dir_target Name of the variable that holds the target directory
MACRO (OPENMS_COPYDIR dir_source dir_target)

  ## deactivated:
  ## (checking if the directory exists does not suffice, because the content (i.e. headers) might have changed due to new patches)
  #if (EXISTS ${${dir_target}})
  # Message(STATUS "Copying ${${dir_source}} --> ${${dir_target}} .. skipped (already exists)")
  #else()
    Message(STATUS "Copying ${${dir_source}} --> ${${dir_target}} .. ")
    exec_program(${CMAKE_COMMAND}
      ARGS -E copy_directory "\"${${dir_source}}\"" "\"${${dir_target}}\""
      OUTPUT_VARIABLE COPYDIR_OUT
      RETURN_VALUE COPY_SUCCESS)

    # logfile
    file(APPEND ${LOGFILE} ${COPYDIR_OUT})

    if (NOT COPY_SUCCESS EQUAL 0)
      message(FATAL_ERROR "Copying ${${dir_source}} --> ${${dir_target}} .. done")
      message(FATAL_ERROR ${COPYDIR_OUT})
    else()
      message(STATUS "Copying ${${dir_source}} --> ${${dir_target}} .. done")
    endif()
  #endif()
ENDMACRO (OPENMS_COPYDIR)

## default header for each library in the log file
## @param libname Name of the lirary that is shown in the library
MACRO (OPENMS_LOGHEADER_LIBRARY libname)
    file(APPEND ${LOGFILE} "\n\r")
    file(APPEND ${LOGFILE} " ------------------------------------------------------------------------\n\r")
    file(APPEND ${LOGFILE} "                              ${libname}\n\r")
    file(APPEND ${LOGFILE} " ------------------------------------------------------------------------\n\r")
    file(APPEND ${LOGFILE} "\n\r")
ENDMACRO (OPENMS_LOGHEADER_LIBRARY)

#############################################################################################
##################################  CLEANING MACROS  ########################################
#############################################################################################

## Removes all installed files from the specified contrib package
## @param libname Name of the package that should be cleaned
MACRO (OPENMS_CLEAN_LIB libname)
  if( FORCE_REBUILD MATCHES "ON") # find & delete all installed
    OPENMS_CLEAN_INSTALLED_LIBS(${libname})
    OPENMS_CLEAN_HEADERS(${libname})
    OPENMS_CLEAN_SOURCE(${libname})
  endif()
ENDMACRO (OPENMS_CLEAN_LIB libname)

# removes all files installed on windows libraries
# @param libname Name of the library that should be removed
MACRO (OPENMS_CLEAN_INSTALLED_LIBS libname)
  string(TOUPPER ${libname} libnameUP)
  # remove libs
  message(STATUS "Removing library files installed by ${libname} .. ")
  if(MSVC)
    file(GLOB_RECURSE LIB_FILES "\"${${libnameUP}_DIR}/*.lib\"" "\"${${libnameUP}_DIR}/*.dll\"")
  else()
    file(GLOB_RECURSE LIB_FILES "\"${${libnameUP}_DIR}/*.a\"" "\"${${libnameUP}_DIR}/*.la\"" "\"${${libnameUP}_DIR}/*.o\"")
  endif()

  set(REMOVED_LIB 0)

  foreach (FTD ${LIB_FILES})
    get_filename_component(RFTD ${FTD} NAME)
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove "\"${CONTRIB_BIN_LIB_DIR}/${RFTD}\""
                    OUTPUT_VARIABLE DELETE_LIB_OUT
                    RESULT_VARIABLE DELETE_LIB_SUCCESS)
    if( NOT DELETE_LIB_SUCCESS EQUAL 0)
      message(STATUS "Removing library files installed by ${libname} .. failed")
      file(APPEND ${LOGFILE} ${DELETE_LIB_OUT})
      message( FATAL_ERROR ${DELETE_LIB_OUT})
    endif()
    set(REMOVED_LIB 1)
  endforeach (FTD ${LIB_FILES})
  if(${REMOVED_LIB})
    message(STATUS "Removing library files installed by ${libname} .. done")
  else()
    message(STATUS "Removing library files installed by ${libname} .. skipped")
  endif()
ENDMACRO (OPENMS_CLEAN_INSTALLED_LIBS libname)

# remove all header files installed by the library
# @param libname Name of the library
MACRO (OPENMS_CLEAN_HEADERS libname)
 string(TOUPPER ${libname} libnameUP)
 # delete header files
  message(STATUS "Removing header files installed by ${libname} .. ")
  set(REMOVED_HEADER 0)
  # handle SVM extra
  if(DEFINED INCLUDE_FILES_${libnameUP})

    foreach(FTD ${INCLUDE_FILES_${libnameUP}})
      if(EXISTS ${FTD})
        execute_process(COMMAND ${CMAKE_COMMAND} -E remove "${FTD}"
                        OUTPUT_VARIABLE DELETE_HEADER_FILE_OUT
                        RESULT_VARIABLE DELETE_HEADER_FILE_SUCCESS)
        if( NOT DELETE_HEADER_FILE_SUCCESS EQUAL 0)
          message(STATUS "Removing header files installed by ${libname} .. failed")
          file(APPEND ${LOGFILE} ${DELETE_SVM_HEADER_OUT})
          message( FATAL_ERROR ${DELETE_SVM_HEADER_OUT})
        endif()
        set(REMOVED_HEADER 1)
      endif()
    endforeach()

  else()
    # remove complete include directory
    if(EXISTS ${INCLUDE_DIR_${libnameUP}})
      execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory "${INCLUDE_DIR_${libnameUP}}"
                      OUTPUT_VARIABLE DELETE_DIR_HEADER_OUT
                      RESULT_VARIABLE DELETE_DIR_HEADER_SUCCESS)
      if( NOT DELETE_DIR_HEADER_SUCCESS EQUAL 0)
        message(STATUS "Removing header files installed by ${libname} .. failed")
        file(APPEND ${LOGFILE} ${DELETE_DIR_HEADER_OUT})
        message( FATAL_ERROR ${DELETE_DIR_HEADER_OUT})
      endif()
      set(REMOVED_HEADER 1)
    endif()
  endif ()

  if(${REMOVED_HEADER})
    message(STATUS "Removing header files installed by ${libname} .. done")
  else()
    message(STATUS "Removing header files installed by ${libname} .. skipped")
  endif()
ENDMACRO(OPENMS_CLEAN_HEADERS libname)

MACRO(OPENMS_CLEAN_SOURCE libname)
  string(TOUPPER ${libname} libnameUP)
  message(STATUS "Removing ${libname} source directory .. ")
  if(EXISTS "${${libnameUP}_DIR}")
    # delete source directory
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory "\"${${libnameUP}_DIR}\""
                    OUTPUT_VARIABLE DELETE_SRC_DIR_OUT
                    RESULT_VARIABLE DELETE_SRC_DIR_SUCCESS)
    if( NOT DELETE_SRC_DIR_SUCCESS EQUAL 0)
      message(STATUS "Removing ${libname} source directory .. failed")
      file(APPEND ${LOGFILE} ${DELETE_SRC_DIR_OUT})
      message( FATAL_ERROR ${DELETE_SRC_DIR_OUT})
    endif()
    message(STATUS "Removing ${libname} source directory .. done")
  else()
    message(STATUS "Removing ${libname} source directory .. skipped")
  endif()

  if(MSVC) #windows: delete intermediate tar
    if(DEFINED ARCHIVE_${libnameUP}_TAR)
      message(STATUS "Removing ${libname} intermediate archives .. ")
      if(EXISTS ${ARCHIVE_${libnameUP}_TAR})
        execute_process(COMMAND ${CMAKE_COMMAND} -E remove "\"${ARCHIVE_${libnameUP}_TAR}\""
                        OUTPUT_VARIABLE DELETE_INTERMEDIATE_TAR_OUT
                        RESULT_VARIABLE DELETE_INTERMEDIATE_TAR_SUCCESS)
        if( NOT DELETE_INTERMEDIATE_TAR_SUCCESS EQUAL 0)
          message(STATUS "Removing ${libname} intermediate archives .. failed")
          file(APPEND ${LOGFILE} ${DELETE_INTERMEDIATE_TAR_OUT})
          message( FATAL_ERROR ${DELETE_INTERMEDIATE_TAR_OUT})
        endif()
        message(STATUS "Removing ${libname} intermediate archives .. done")
      else()
        message(STATUS "Removing ${libname} intermediate archives .. skipped")
      endif()
    endif()
  endif(MSVC)
ENDMACRO(OPENMS_CLEAN_SOURCE libname)

#############################################################################################
####################################  COPY MACROS  ##########################################
#############################################################################################


## copies all created lib files from a specific library to the
## contrib/lib dir -- is only active under MSVC
## @param libname Name of the library
MACRO(OPENMS_COPY_LIBS libname)
  ## TODO: Revise this and look for another solution
  ## copy *.lib, *.dll  etc to ./lib
  string(TOUPPER ${libname} libnameUP)
  if (MSVC)
    message(STATUS "Copying ${libname} libraries .. ")
    file(GLOB_RECURSE LIB_FILES "${${libnameUP}_DIR}/*.lib" "${${libnameUP}_DIR}/*.dll")

    foreach (FTC ${LIB_FILES})
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${FTC}" "${CONTRIB_BIN_LIB_DIR}"
                      OUTPUT_VARIABLE COPY_LIB_OUT
                      RESULT_VARIABLE COPY_LIB_SUCCESS)
      if( NOT COPY_LIB_SUCCESS EQUAL 0)
        message( STATUS "Copying ${libname} libraries (${FTC}) .. failed")
        file(APPEND ${LOGFILE} ${COPY_LIB_OUT})
        message( FATAL_ERROR ${COPY_LIB_OUT})
      endif()
    endforeach()
    message(STATUS "Copying ${libname} libraries .. done")
  endif()
ENDMACRO(OPENMS_COPY_LIBS libname)

########
########
