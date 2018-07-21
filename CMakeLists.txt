cmake_minimum_required(VERSION 2.8)
project(libqtree)

option(ENABLE_INSTRUMENTATION "" OFF)
if(ENABLE_INSTRUMENTATION)
  add_definitions(-DHAS_INSTR)
endif()

find_package(GTest REQUIRED)
find_package(GMock REQUIRED)

enable_testing()

add_library(qdtree src/utils.cxx)
target_include_directories(qdtree
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(unit-tests
  test/octree.cxx
  test/qdtree-test.cxx)

target_include_directories(unit-tests
  PRIVATE {CMAKE_CURRENT_SOURCE_DIR}/test
          ${GTEST_INCLUDE_DIRS}
          ${GMOCK_INCLUDE_DIRS})

target_link_libraries(unit-tests
  qdtree
  GTest::GTest
  ${GMOCK_BOTH_LIBRARIES})

gtest_add_tests(unit-tests "" AUTO)
