# =========================================================
#
# Some cmake macros to deal with tests.
# =========================================================
MACRO(BEGIN_TEST _D)
  SET(_CURRENT_TEST_DIRECTORY ${_D})
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_D})

  # find and copy data files : *.mat, *.dat and *.xml, and etc.
  FILE(GLOB_RECURSE _DATA_FILES 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_D}
    *.mat
    *.dat
    *.hdf5
    *.xml
    *.DAT
    *.INI)
  FOREACH(_F ${_DATA_FILES})
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${_D}/${_F}
      ${CMAKE_CURRENT_BINARY_DIR}/${_D}/${_F} COPYONLY)
  ENDFOREACH(_F ${_DATA_FILES})

  IF(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_D}-utils)
    FILE(GLOB_RECURSE TEST_UTILS_SOURCES_TMP ${CMAKE_CURRENT_SOURCE_DIR}/${_D}-utils/*.c)
    set(TEST_UTILS_SOURCES ${TEST_UTILS_SOURCES} ${TEST_UTILS_SOURCES_TMP})
    SET(${_CURRENT_TEST_DIRECTORY}_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/${_D}-utils)
  ELSE()
    SET(${_CURRENT_TEST_DIRECTORY}_INCLUDE_DIRECTORIES)
  ENDIF()

  # configure test CMakeLists.txt (needed for a chdir before running test)
  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/CMakeListsForTests.cmake 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/CMakeLists.txt @ONLY)

  SET(_EXE_LIST_${_CURRENT_TEST_DIRECTORY})
ENDMACRO(BEGIN_TEST _D)

# Tests
MACRO(BEGIN_TEST2 _D)
  SET(_CURRENT_TEST_DIRECTORY ${_D})
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_D})

  # find and copy data files : *.mat, *.dat and *.xml, and etc.
  FILE(GLOB_RECURSE _DATA_FILES 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_D}
    *.mat 
    *.dat
    *.hdf5
    *.xml
    *.DAT
    *.INI)

  IF(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_D}-utils)
    FILE(GLOB_RECURSE TEST_UTILS_SOURCES_TMP ${CMAKE_CURRENT_SOURCE_DIR}/${_D}-utils/*.[ch])
    set(TEST_UTILS_SOURCES ${TEST_UTILS_SOURCES} ${TEST_UTILS_SOURCES_TMP})
    SET(${_CURRENT_TEST_DIRECTORY}_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/${_D}-utils)
  ELSE()
    SET(${_CURRENT_TEST_DIRECTORY}_INCLUDE_DIRECTORIES)
  ENDIF()

  FOREACH(_F ${_DATA_FILES})
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${_D}/${_F} ${CMAKE_CURRENT_BINARY_DIR}/${_D}/${_F} COPYONLY)
  ENDFOREACH(_F ${_DATA_FILES})

  # configure test CMakeLists.txt (needed for a chdir before running test)
  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/CMakeListsForTestsv2.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/CMakeLists.txt @ONLY)

  SET(_EXE_LIST_${_CURRENT_TEST_DIRECTORY})

ENDMACRO(BEGIN_TEST2 _D _L)

# Declaration of a siconos test
MACRO(NEW_TEST)
  CAR(_EXE ${ARGV})
  CDR(_SOURCES ${ARGV})
  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${_EXE})
  SET(${_EXE}_FSOURCES)
  FOREACH(_F ${_SOURCES})
    LIST(APPEND ${_EXE}_FSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${_CURRENT_TEST_DIRECTORY}/${_F})
  ENDFOREACH(_F ${_SOURCES})
 
  IF(TEST_MAIN)
    LIST(APPEND ${_EXE}_FSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_MAIN})
  ENDIF(TEST_MAIN)

ENDMACRO(NEW_TEST)

# Removal of a siconos test (test fails or takes forever)
# Warning: the test is still compiled
MACRO(RM_TEST)
  CAR(_EXE ${ARGV})
  CDR(_SOURCES ${ARGV})
  LIST(REMOVE_ITEM _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${_EXE})
ENDMACRO(RM_TEST)

MACRO(RM_TEST2)
  SET(TEST_SOLVER ${ARGV0})
  SET(TEST_DATA ${ARGV1})

  SET(TEST_SBM ${ARGV2})
  SET(TEST_SBM_C "_SBM")
  IF(NOT DEFINED TEST_SBM)
    SET(TEST_SBM 0)
    SET(TEST_SBM_C "")
  ENDIF(NOT DEFINED TEST_SBM)

  STRING(REGEX REPLACE SICONOS_ "" TEST_SOLVER_NAME ${TEST_SOLVER})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})
  SET(TEST_EXE ${TEST_SOLVER_NAME}${TEST_SBM_C})
  SET(TEST_NAME "test-${TEST_SOLVER_NAME}-${TEST_DATA_NAME}${TEST_SBM_C}")

  LIST(REMOVE_ITEM ${TEST_EXE}_DATA_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_DATA})
  LIST(REMOVE_ITEM ${TEST_EXE}_NAME_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})
ENDMACRO(RM_TEST2)

MACRO(NEW_FC_TEST)
  SET(TEST_SOLVER ${ARGV0})
  SET(TEST_DATA ${ARGV1})
  
  SET(TEST_TOLERANCE ${ARGV2})
  IF(NOT DEFINED TEST_TOLERANCE)
    SET(TEST_TOLERANCE 0)
  ENDIF(NOT DEFINED TEST_TOLERANCE)
  
  SET(TEST_MAXITER ${ARGV3})
  IF(NOT DEFINED TEST_MAXITER)
    SET(TEST_MAXITER 0)
  ENDIF(NOT DEFINED TEST_MAXITER)
  
  SET(TEST_INTERNAL_SOLVER ${ARGV4})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER)
    SET(TEST_INTERNAL_SOLVER 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER)
  
  SET(TEST_INTERNAL_SOLVER_TOLERANCE ${ARGV5})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER_TOLERANCE)
    SET(TEST_INTERNAL_SOLVER_TOLERANCE 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER_TOLERANCE)
  
  SET(TEST_INTERNAL_SOLVER_MAXITER ${ARGV6})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER_MAXITER)
    SET(TEST_INTERNAL_SOLVER_MAXITER 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER_MAXITER)
  
  STRING(REGEX REPLACE SICONOS_FRICTION_ "" TEST_SOLVER_NAME ${TEST_SOLVER})
  STRING(REGEX REPLACE SICONOS_FRICTION_ "" TEST_INTERNAL_SOLVER_NAME ${TEST_INTERNAL_SOLVER})
  STRING(REGEX REPLACE "0" "" TEST_INTERNAL_SOLVER_NAME ${TEST_INTERNAL_SOLVER_NAME})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})

  SET(TEST_NAME "test-${TEST_SOLVER_NAME}${TEST_INTERNAL_SOLVER_NAME}-${TEST_DATA_NAME}")

  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/fctest.c.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

  SET(${TEST_NAME}_FSOURCES)

  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})
  LIST(APPEND ${TEST_NAME}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

ENDMACRO(NEW_FC_TEST)

MACRO(NEW_PB_TEST)
  SET(FILE_TO_CONF ${ARGV0})
  SET(TEST_SOLVER ${ARGV1})
  SET(TEST_DATA ${ARGV2})

  SET(TEST_SBM ${ARGV3})
  SET(TEST_SBM_C "_SBM")
  IF(NOT DEFINED TEST_SBM)
    SET(TEST_SBM 0)
    SET(TEST_SBM_C "")
  ENDIF(NOT DEFINED TEST_SBM)


  STRING(REGEX REPLACE SICONOS_ "" TEST_SOLVER_NAME ${TEST_SOLVER})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})

  SET(TEST_EXE ${TEST_SOLVER_NAME}${TEST_SBM_C})
  SET(TEST_NAME "test-${TEST_SOLVER_NAME}-${TEST_DATA_NAME}${TEST_SBM_C}")


  LIST(FIND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_EXE} ALREADY_CONF)
  IF(ALREADY_CONF EQUAL -1)
    SET(${TEST_EXE}_FSOURCES)

    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/${FILE_TO_CONF}${TEST_SBM_C}.c.in
      ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_EXE}.c)
    LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_EXE})
    LIST(APPEND ${TEST_EXE}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_EXE}.c)
    SET(${TEST_EXE}_DATA_LIST_${_CURRENT_TEST_DIRECTORY})
    SET(${TEST_EXE}_NAME_LIST_${_CURRENT_TEST_DIRECTORY})
  ENDIF(ALREADY_CONF EQUAL -1)

  LIST(APPEND ${TEST_EXE}_DATA_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_DATA})
  LIST(APPEND ${TEST_EXE}_NAME_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})

ENDMACRO(NEW_PB_TEST)

MACRO(NEW_LCP_TEST)
  NEW_PB_TEST(lcptest ${ARGV})
ENDMACRO(NEW_LCP_TEST)

MACRO(NEW_RELAY_TEST)
  NEW_PB_TEST(relaytest ${ARGV})
ENDMACRO(NEW_RELAY_TEST)

MACRO(NEW_NCP_TEST)
  SET(FILE_TO_CONF ${ARGV0})
  SET(TEST_SOLVER ${ARGV1})
  STRING(REGEX REPLACE SICONOS_ "" TEST_SOLVER_NAME ${TEST_SOLVER})
  SET(TEST_EXE ${TEST_SOLVER_NAME}-${FILE_TO_CONF})
  SET(TEST_NAME "test-${TEST_SOLVER_NAME}-${FILE_TO_CONF}")

  SET(${TEST_EXE}_FSOURCES)

  SET(SOLVER_ID ${TEST_SOLVER})
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${_CURRENT_TEST_DIRECTORY}/${FILE_TO_CONF}.c.in
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_EXE}.c)
  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_EXE})
  LIST(APPEND ${TEST_EXE}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_EXE}.c)
  SET(${TEST_EXE}_DATA_LIST_${_CURRENT_TEST_DIRECTORY} )
  SET(${TEST_EXE}_NAME_LIST_${_CURRENT_TEST_DIRECTORY})
ENDMACRO(NEW_NCP_TEST)

MACRO(NEW_LS_TEST)
  SET(TEST_SOLVER ${ARGV0})
  SET(TEST_DATA ${ARGV1})

  SET(TEST_SBM ${ARGV2})
  SET(TEST_SBM_C "SBM")
  IF(NOT DEFINED TEST_SBM)
    SET(TEST_SBM 0)
    SET(TEST_SBM_C "")
  ENDIF(NOT DEFINED TEST_SBM)
    
  STRING(REGEX REPLACE SICONOS_ "" TEST_SOLVER_NAME ${TEST_SOLVER})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})

  SET(TEST_NAME "test-${TEST_SOLVER_NAME}-${TEST_DATA_NAME}${TEST_SBM_C}")

  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/lstest.c.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

  SET(${TEST_NAME}_FSOURCES)

  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})
  LIST(APPEND ${TEST_NAME}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

ENDMACRO(NEW_LS_TEST)

MACRO(NEW_GMP_TEST)
  SET(TEST_SOLVER ${ARGV0})
  SET(TEST_DATA ${ARGV1})
  SET(TEST_GMP_REDUCED 1)
  
  SET(TEST_TOLERANCE ${ARGV2})
  IF(NOT DEFINED TEST_TOLERANCE)
    SET(TEST_TOLERANCE 0)
  ENDIF(NOT DEFINED TEST_TOLERANCE)
  
  SET(TEST_MAXITER ${ARGV3})
  IF(NOT DEFINED TEST_MAXITER)
    SET(TEST_MAXITER 0)
  ENDIF(NOT DEFINED TEST_MAXITER)
  
  SET(TEST_INTERNAL_SOLVER ${ARGV4})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER)
    SET(TEST_INTERNAL_SOLVER 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER)
  
  SET(TEST_INTERNAL_SOLVER_TOLERANCE ${ARGV5})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER_TOLERANCE)
    SET(TEST_INTERNAL_SOLVER_TOLERANCE 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER_TOLERANCE)
  
  SET(TEST_INTERNAL_SOLVER_MAXITER ${ARGV6})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER_MAXITER)
    SET(TEST_INTERNAL_SOLVER_MAXITER 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER_MAXITER)

  SET(TEST_GMP_REDUCED ${ARGV7})
  IF(NOT DEFINED TEST_GMP_REDUCED)
    SET(TEST_GMP_REDUCED 0)
  ENDIF(NOT DEFINED TEST_GMP_REDUCED)
  
  STRING(REGEX REPLACE SICONOS_FRICTION_ "" TEST_SOLVER_NAME "REDUCED" ${TEST_GMP_REDUCED}_ ${TEST_SOLVER})
  STRING(REGEX REPLACE SICONOS_FRICTION_ "" TEST_INTERNAL_SOLVER_NAME ${TEST_INTERNAL_SOLVER})
  STRING(REGEX REPLACE "0" "" TEST_INTERNAL_SOLVER_NAME ${TEST_INTERNAL_SOLVER_NAME})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})

  SET(TEST_NAME "test-GMP-${TEST_SOLVER_NAME}${TEST_INTERNAL_SOLVER_NAME}-${TEST_DATA_NAME}")

  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/gmptest.c.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

  SET(${TEST_NAME}_FSOURCES)

  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})
  LIST(APPEND ${TEST_NAME}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

ENDMACRO(NEW_GMP_TEST)

# add subdirs (i.e. CMakeLists.txt generated for tests) to the build
MACRO(END_TEST)
  ADD_SUBDIRECTORY(${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY} ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY})
ENDMACRO(END_TEST)


# Build plugins required for python tests
macro(build_plugin plug)
  get_filename_component(plug_name ${plug} NAME_WE)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests/plugins/)
  add_library(${plug_name} MODULE ${plug})
  set_property(TARGET ${plug_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/tests)
  set_target_properties(${plug_name} PROPERTIES PREFIX "")
  add_dependencies(${COMPONENT} ${plug_name})
endmacro()

# ----------------------------------------
# Prepare python tests for the current
# component
# Usage:
#   build_python_tests(path_to_tests)
#
# path_to_tests is relative to the current source dir.
# Most of the time, path_to_tests = 'tests'.
# For instance, in mechanics, tests are called in CMakeLists.txt
# in swig, current source dir is thus source_dir/mechanics/swig
# and source_dir/mechanics/swig/tests contains all the python files for tests.
#
# This routine copy the directory of tests to binary dir to allow 'py.test' run in the build.
# 
# binary dir will then look like :
# wrap/siconos
# wrap/siconos/mechanics
# wrap/siconos/mechanics/tests
#
# and running py.tests in wrap dir will end up with a run of
# all mechanics tests.
macro(build_python_tests)
  if(WITH_${COMPONENT}_TESTING)
    # copy data files
    #file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${_D}/data ${SICONOS_SWIG_ROOT_DIR}/${_D}/data)
    
    # build plugins, if any
    # Note : all built libraries are saved in SICONOS_SWIG_ROOT_DIR/plugins
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/plugins)
      file(GLOB plugfiles ${CMAKE_CURRENT_SOURCE_DIR}/tests/plugins/*.cpp)
      foreach(plug ${plugfiles})
	build_plugin(${plug})
      endforeach()
    endif()
    
    # copy test dir to binary dir (inside siconos package)
    # ---> allows py.test run in binary dir
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)
      file(GLOB data4tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/tests/data
	${CMAKE_CURRENT_SOURCE_DIR}/tests/data/*)
      foreach(datafile ${data4tests})
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/data/${datafile}
	  ${SICONOS_SWIG_ROOT_DIR}/tests/data/${datafile} COPYONLY)
      endforeach()
    endif()
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/CAD)
      file(GLOB data4tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/tests/CAD
	${CMAKE_CURRENT_SOURCE_DIR}/tests/CAD/*)
      foreach(datafile ${data4tests})
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/CAD/${datafile}
	  ${SICONOS_SWIG_ROOT_DIR}/tests/CAD/${datafile} COPYONLY)
      endforeach()
    endif()
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests)
      file(GLOB testfiles ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.py)
      foreach(excluded_test ${${COMPONENT}_python_excluded_tests})
	list(REMOVE_ITEM testfiles ${excluded_test})
      endforeach()
      foreach(file ${testfiles})
	get_filename_component(testname ${file} NAME_WE)
	get_filename_component(exename ${file} NAME)
	# Each file is copy to siconos/tests.
	# Maybe we can create a 'tests' dir for each subpackage?
	# --> Easier to deal with plugins and data if only one package
	configure_file(${file} ${SICONOS_SWIG_ROOT_DIR}/tests COPYONLY)
	set(name "python_${testname}")
	set(exename ${SICONOS_SWIG_ROOT_DIR}/tests/${exename})
	add_python_test(${name}, ${exename})
      endforeach()
    endif()
  endif()
endmacro()

# Declaration of a siconos test based on python bindings
macro(add_python_test test_name test_file)
  add_test(${test_name} ${PYTHON_EXECUTABLE} ${TESTS_RUNNER} "${pytest_opt}" ${DRIVE_LETTER}${test_file})
  #    WORKING_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/tests)
  set_tests_properties(${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "FAILURE;Exception;[^x]failed;ERROR;Assertion")
  set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=$ENV{PYTHONPATH}:${CMAKE_BINARY_DIR}/wrap")
  set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/wrap/siconos/tests") # for plugins
endmacro()
