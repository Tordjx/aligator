#
# Copyright (C) 2022-2023 LAAS-CNRS, INRIA
#

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

macro(ADD_TEST_CFLAGS test_name flag)
  set_property(
    TARGET ${test_name}
    APPEND_STRING
    PROPERTY COMPILE_FLAGS " ${flag}"
  )
endmacro()

file(GLOB TEST_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
message(STATUS "Test headers: ${TEST_HEADERS}")

function(_add_test_prototype name prefix dependencies)
  set(test_name "${PROJECT_NAME}-test-cpp-${prefix}${name}")
  set(test_file ${name}.cpp)
  message(STATUS "Adding cpp test: ${test_file} (${test_name})")

  ADD_UNIT_TEST(${test_name} ${test_file} ${TEST_HEADERS})
  set_target_properties(${test_name} PROPERTIES LINKER_LANGUAGE CXX)
  set_standard_output_directory(${test_name})
  target_include_directories(${test_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  # define macros required by boost_test
  # see: https://www.boost.org/doc/libs/1_78_0/libs/test/doc/html/boost_test/usage_variants.html
  # define module name, replace '-' by '_'
  set(MODULE_NAME "${name}Test")
  string(REGEX REPLACE "-" "_" MODULE_NAME ${MODULE_NAME})

  ADD_TEST_CFLAGS(${test_name} "-DBOOST_TEST_DYN_LINK")
  ADD_TEST_CFLAGS(${test_name} "-DBOOST_TEST_MODULE=${MODULE_NAME}")

  target_link_libraries(${test_name} PUBLIC ${dependencies})
  target_link_libraries(${test_name} PRIVATE Boost::unit_test_framework)
endfunction()

function(add_aligator_test name)
  _add_test_prototype(${name} "" ${PROJECT_NAME})
endfunction()

set(
  TEST_NAMES
  continuous
  costs
  integrators
  lqr
  problem
  solver-storage
  utils
)

if(PROXSUITE_NLP_WITH_PINOCCHIO_SUPPORT)
  list(APPEND TEST_NAMES cycling mpc-cycle)
endif()

if(BUILD_WITH_PINOCCHIO_SUPPORT AND PINOCCHIO_V3)
  list(APPEND TEST_NAMES forces)
endif()

foreach(test_name ${TEST_NAMES})
  add_aligator_test(${test_name})
endforeach()

if(CHECK_RUNTIME_MALLOC)
  add_aligator_test(nomalloc)
endif()

add_subdirectory(gar)

# PYTHON TESTS
if(BUILD_PYTHON_INTERFACE)
  add_subdirectory(python)
endif()

if(BUILD_CROCODDYL_COMPAT)
  add_subdirectory(compat)
endif()
