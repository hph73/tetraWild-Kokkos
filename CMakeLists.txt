cmake_minimum_required (VERSION 3.10)
project (Testing)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(KokkosKernels REQUIRED)
include(./common.cmake)
find_package(GMPfTetWild)
IF(NOT ${GMP_FOUND})
	MESSAGE(FATAL_ERROR "Cannot find GMP")
ENDIF()

include_directories(
  ./include
  ./include/eigen
  ./include/libigl
)

set(
  SOURCES
  ./src/main.cpp
  ./src/tetraWildFunctions.cpp
)

add_executable(main ${SOURCES})
target_link_libraries(
  main
  ${CMAKE_SOURCE_DIR}/lib/libFloatTetwild.a
  ${CMAKE_SOURCE_DIR}/lib/libgeogram.a
  ${CMAKE_SOURCE_DIR}/lib/libfmt.a
  ${GMP_LIBRARIES}
  Kokkos::kokkoskernels
)
