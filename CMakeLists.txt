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

add_library(quadtree STATIC test/quadtree.cxx)
target_include_directories(quadtree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test)
target_link_libraries(quadtree PUBLIC qdtree)

add_executable(unit-tests test/qdtree-test.cxx)
target_link_libraries(unit-tests
  PRIVATE quadtree ${GMOCK_MAIN_LIBRARIES} Threads::Threads)
gtest_add_tests(unit-tests "" AUTO)

find_package(Boost COMPONENTS program_options REQUIRED)

add_executable(benchmark test/benchmark.cxx)
target_include_directories(benchmark PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(benchmark PRIVATE quadtree ${Boost_LIBRARIES})
