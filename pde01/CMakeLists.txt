# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project(pde01 VERSION 1.0 DESCRIPTION "simple tests for vs to cmake conversion")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(pde01 
	doubleheateq.cpp
	ScopedWindowsPerformanceTimer.cpp
	stdafx.cpp
	testQuantlib.cpp
)

target_link_libraries(pde01
	Approximation 
	libs 
	PointReductionLib_use_dll 	
)

set(UTIL_FILES 
		"${CMAKE_CURRENT_SOURCE_DIR}/QuadraturePoints12.txt" 
		"${CMAKE_CURRENT_SOURCE_DIR}/QuadraturePoints14.txt" 
		"${CMAKE_CURRENT_SOURCE_DIR}/QuadraturePoints16.txt" 
		"${CMAKE_CURRENT_SOURCE_DIR}/QuadraturePoints18.txt" 
		"${CMAKE_CURRENT_SOURCE_DIR}/QuadraturePoints10.txt" 
		"${CMAKE_CURRENT_SOURCE_DIR}/QuadraturePoints20.txt" 
)

# Set the output directory for the executable
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Copy the library files to the executable directory
#file(COPY ${LIB_FILES} DESTINATION "${EXECUTABLE_OUTPUT_PATH}")
#file(COPY ${UTIL_FILES} DESTINATION${EXECUTABLE_OUTPUT_PATH})

add_custom_command(TARGET pde01 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${UTIL_FILES}
        $<TARGET_FILE_DIR:pde01>)

add_custom_command(TARGET pde01 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${INTEL_OMP_DLL}
        $<TARGET_FILE_DIR:pde01>)

target_include_directories(pde01 PRIVATE
	${CMAKE_SOURCE_DIR}/Approximation
	${CMAKE_SOURCE_DIR}/libs
	${CMAKE_SOURCE_DIR}/PointReductionLib_use_dll
	${CMAKE_SOURCE_DIR}/recombine
	${CMAKE_SOURCE_DIR}/CommandLine
	${CMAKE_SOURCE_DIR}/Morton
	${BOOST_INCLUDE_DIRS}
	${QUANTLIB_INCLUDE_DIRS}
)

target_link_directories(pde01 PRIVATE
	${BOOST_LIBRARIES}
	${QUANTLIB_LIBRARIES}
	${MKL_LIBRARIES}
)

target_link_libraries(pde01
	Approximation
	libs
	PointReductionLib_use_dll
	recombine
	CommandLine
	Morton
)

