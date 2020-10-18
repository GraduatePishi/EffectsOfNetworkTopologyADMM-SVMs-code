Required Tools/Software
1.CMake
2.MPICH


TO DO:
When cloning the code for the first time from git, 
in Windows:
in order Gtest(GoogleTest) works, "configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)" should be active,
and it should exists a CMAKELists.txt.in file in the folder with the content of 

cmake_minimum_required(VERSION 2.8.2)
project(googletest-download NONE)
include(ExternalProject)
ExternalProject_Add(googletest
  GIT_REPOSITORY    https://github.com/google/googletest.git
  GIT_TAG           master
  SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-src"
  BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)
