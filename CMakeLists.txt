# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_GENERATOR_TOOLSET "v120") # choose msvc 2013

project ("top")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#message(STATUS "Compiler: ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_VERSION}")


set(MKL_INCLUDE_DIRS "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include")
set(MKL_LIBRARIES  "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64_win")
set(BOOST_INCLUDE_DIRS "C:/Program Files/boost/boost_1_66_0")
set(BOOST_LIBRARIES "C:/Program Files/boost/boost_1_66_0/stage_x64/lib" )
set(QUANTLIB_INCLUDE_DIRS "C:/Users/pan_w/Downloads/QuantLib-1.21")
set(QUANTLIB_LIBRARIES "C:/Users/pan_w/Downloads/QuantLib-1.21/lib")
set(INTEL_OMP_DLL "${CMAKE_CURRENT_SOURCE_DIR}/MKL64/libiomp5md.dll")
set(INTEL_OMP_LIB "${CMAKE_CURRENT_SOURCE_DIR}/MKL64/libiomp5md.lib")

function(check_exists filepath)
   if(NOT EXISTS "${filepath}")
      message(FATAL_ERROR "${filepath} does not exist")
   endif()
endfunction()

# Include sub-projects.
add_subdirectory(Morton)		# static lib - builds
add_subdirectory(Approximation) # static lib - builds
add_subdirectory(CommandLine)	# static lib - builds
add_subdirectory(libs)			# static lib - builds
add_subdirectory(PointReductionLib_use_dll) # static lib - builds
add_subdirectory(recombine)		# static lib - builds
add_subdirectory(pde01)			#executable


add_dependencies(Approximation Morton)
add_dependencies(libs CommandLine)
add_dependencies(PointReductionLib_use_dll Morton recombine libs)
add_dependencies(pde01 Approximation libs PointReductionLib_use_dll)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT pde01)
