# CMakeList.txt : CMake project for MachineLearningPackage_CMAKE, include source and define
# project specific logic here.
cmake_minimum_required (VERSION 3.11)
project(MachineLearningExec)
#################################
########## Google Test ########## I have no changes in this part
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
if(WIN32)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
		RESULT_VARIABLE result
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
	if(result)
		 message(FATAL_ERROR "CMake step for googletest failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
		RESULT_VARIABLE result
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
	if(result)
		message(FATAL_ERROR "Build step for googletest failed: ${result}")
	endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
	add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
					${CMAKE_CURRENT_BINARY_DIR}/googletest-build
					EXCLUDE_FROM_ALL)
include_directories("${CMAKE_CURRENT_BINARY_DIR}/googletest-src/googlemock/include")
# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
	if (CMAKE_VERSION VERSION_LESS 2.8.11)
		include_directories("${gtest_SOURCE_DIR}/include")
	endif()
endif()
############ End Of Google Test #####################
#build properties
set(CMAKE_CXX_STANDARD 17)

#Preparing MPI requirements
find_package(MPI REQUIRED)

if(MPI_FOUND) #aparently ENABLE_MPI does not enable the MPI in windows
	message("*******  MPI is ON ******* ")
	set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} "-O2") 
	set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}) 
	include_directories(SYSTEM ${MPI_INCLUDE_PATH}) 
else()
	message("!!!!!!!!!!!  MPI NOT FOUND !!!!!!!!!!!!")
endif(MPI_FOUND)

#addreses to libraries and headerfiles
add_subdirectory(external/fmt)


find_library(NLOPTPATH NAMES nlopt HINTS ENV LIBRARY_PATH)
if(NOT NLOPTPATH AND WIN32)
	message("Not finds NLOPT and builds it in windows!!!!!!")
	find_library(NLOPTPATH libnlopt-0 external/nlopt/) 
endif()
if(UNIX)
	message("*******  Builds as linux ******* ")
	set(CMAKE_BUILD_TYPE Relase)
	link_libraries(stdc++fs)
	SET(CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0")
	find_library(STDCPPFSLIB REQUIRE stdc++fs HINTS ENV LIBRARY_PATH)
	include_directories(/apps/Hebbe/software/Core/icc/2019.1.144-GCC-8.2.0-2.31.1/compilers_and_libraries_2019.1.144/linux/pstl/include) 
endif()
message("Flags")
message(${CMAKE_CXX_FLAGS_RELEASE})
include_directories(external)

#build object files
file(GLOB MACHINELEARNING_FILES internal/ADMM_SVM/src/*.h internal/ADMM_SVM/src/*.cpp)
add_executable(MachineLearningExec ${MACHINELEARNING_FILES})
if(WIN32)
	target_link_libraries(MachineLearningExec gtest_main)
	add_test(NAME TestFunctions COMMAND MachineLearningExec)
endif()

#link
if(WIN32)
	message("******* Builds MachineLearningExec in Windows ******* ")
	target_link_libraries(MachineLearningExec fmt ${NLOPTPATH} ${STDCPPFSLIB} ${MPI_LIBRARIES})
else()
	message("******* Builds MachineLearningExec in Linux ******* (2)")
	if (MPI_FOUND)
	message("******* MPI is ON in Linux******* ")
		target_link_libraries(MachineLearningExec fmt ${NLOPTPATH} ${MPI_LIBRARIES})
	endif()
endif()

#post process
if(WIN32)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/external/nlopt/libnlopt-0.dll ${CMAKE_CURRENT_BINARY_DIR}/libnlopt-0.dll COPYONLY)
endif()
install(TARGETS MachineLearningExec DESTINATION bin)
message(STATUS "Build directory "  ${CMAKE_CURRENT_BINARY_DIR})
