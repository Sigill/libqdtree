cmake_minimum_required(VERSION 2.8)
project(libqtree)

option(ENABLE_INSTRUMENTATION "Enable internal log messages" OFF)
if(ENABLE_INSTRUMENTATION)
  add_definitions(-DHAS_INSTR)
endif()

option(BUILD_TEST "Build unit tests and benchmarks" OFF)

add_library(qdtree src/utils.cxx)
target_include_directories(qdtree
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(BUILD_TEST)
  find_package(GTest REQUIRED)
  find_package(GMock REQUIRED)
  enable_testing()

  add_library(quadtree STATIC test/quadtree.cxx)
  target_include_directories(quadtree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test)
  target_link_libraries(quadtree PUBLIC qdtree)

  add_executable(unit-tests test/qdtree-test.cxx)
  target_link_libraries(unit-tests
    PRIVATE quadtree ${GMOCK_MAIN_LIBRARIES} Threads::Threads)
  gtest_add_tests(unit-tests "" AUTO)

  find_package(Boost COMPONENTS program_options REQUIRED)

  include(ExternalProject)

  set(RANG_ROOT ${CMAKE_BINARY_DIR}/thirdparty/rang)

  ExternalProject_Add(rang_external
                      PREFIX ${RANG_ROOT}
                      GIT_REPOSITORY "https://github.com/agauniyal/rang.git"
                      GIT_TAG "49505595c6941ad58788ade0ecdff26aeb5b7a9a"
                      CONFIGURE_COMMAND ""
                      BUILD_COMMAND ""
                      INSTALL_COMMAND "")

  ExternalProject_Get_Property(rang_external SOURCE_DIR)
  set(RANG_INCLUDE_DIR ${SOURCE_DIR}/include)
  file(MAKE_DIRECTORY ${RANG_INCLUDE_DIR})

  add_library(rang INTERFACE IMPORTED)
  set_target_properties(rang PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${RANG_INCLUDE_DIR})
  add_dependencies(rang rang_external)

  add_executable(benchmark test/benchmark.cxx)
  target_include_directories(benchmark
    PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(benchmark
    PRIVATE quadtree ${Boost_LIBRARIES} rang)
endif()
