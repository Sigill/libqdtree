cmake_minimum_required(VERSION 2.8)
project(QDTree)

set(QDTREE_VERSION "0.1.0")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(ENABLE_INSTRUMENTATION "Enable internal log messages" OFF)
if(ENABLE_INSTRUMENTATION)
  add_definitions(-DHAS_INSTR)
endif()

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

set(QDTREE_PUBLIC_HEADERS
  include/qdtree/infix_iterator.hxx
  include/qdtree/node_base.hxx
  include/qdtree/node_base.h
  include/qdtree/singlenode.hxx
  include/qdtree/singlenode.h
  include/qdtree/listnode.hxx
  include/qdtree/listnode.h
  include/qdtree/qdtree.h
  include/qdtree/qdtree.hxx
  include/qdtree/utils.h
  include/qdtree/utils.hxx
  )

add_library(qdtree src/utils.cxx)
target_include_directories(qdtree
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>)

# Define public headers to get them automatically installed.
set_target_properties(qdtree PROPERTIES
  PUBLIC_HEADER "${QDTREE_PUBLIC_HEADERS}")

set(QDTREE_EXPORT QDTreeTargets)

# Install the lib and its headers. Flag it for export.
install(TARGETS qdtree
  EXPORT ${QDTREE_EXPORT}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qdtree
  )

# Create the export file for the build tree.
export(TARGETS qdtree FILE "${PROJECT_BINARY_DIR}/${QDTREE_EXPORT}.cmake")
# Create the export file for the install tree.
install(EXPORT ${QDTREE_EXPORT}
  DESTINATION "${CMAKE_INSTALL_PREFIX}")

# The config files could be generated using a macro (see [1]), but it seems
# to be more complex than to do it manually (see [2]).
# [1] https://cmake.org/cmake/help/latest/module/CMakePackageConfigHelpers.html
# [2] https://gitlab.kitware.com/cmake/community/wikis/doc/tutorials/How-to-create-a-ProjectConfig.cmake-file

# Create the config file for the build tree.
set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
configure_file(cmake/QDTreeConfig.cmake.in
  "${PROJECT_BINARY_DIR}/QDTreeConfig.cmake" @ONLY)

# Create the config file for the install tree.
set(CONF_INCLUDE_DIRS "\${PACKAGE_CMAKE_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
configure_file(cmake/QDTreeConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/QDTreeConfig.cmake" @ONLY)

# Create the package version file.
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/QDTreeConfigVersion.cmake"
  VERSION "${QDTREE_VERSION}"
  COMPATIBILITY SameMajorVersion)

# Install both the config and version file.
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/QDTreeConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/QDTreeConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_PREFIX}")

include(CTest)

if(BUILD_TESTING)
  find_package(GTest REQUIRED NO_MODULE) # Use NO_MODULE to use GTestConfig.cmake instead of FindGTest.cmake
  find_package(Boost COMPONENTS program_options REQUIRED)
  find_package(foonathan_memory REQUIRED)
  find_package(benchmark REQUIRED)


  enable_testing()


  add_library(quadtree STATIC test/quadtree.cxx)
  target_include_directories(quadtree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test)
  target_link_libraries(quadtree PRIVATE qdtree)

  add_executable(unit-tests test/unit-tests.cxx)
  target_link_libraries(unit-tests
    PRIVATE quadtree
            GTest::gmock_main
            Threads::Threads
            foonathan_memory)
  add_test(unit-tests unit-tests)

  add_executable(benchmark test/benchmark.cxx)
  target_link_libraries(benchmark
    PRIVATE quadtree ${Boost_LIBRARIES} benchmark::benchmark_main)

  add_executable(benchmark2 test/benchmark2.cxx)
  target_include_directories(benchmark2 PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(benchmark2
    PRIVATE qdtree ${Boost_LIBRARIES} foonathan_memory benchmark::benchmark_main)
endif()
