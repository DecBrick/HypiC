cmake_minimum_required(VERSION 3.20)

project(HypiC VERSION 0.0.0 LANGUAGES CXX)


#===============================================================================
# Basic CMake configuration
#===============================================================================

# Disable in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR
    "In-source builds are not allowed."
    " Create a separate directory for build files and delete CMakeCache.txt.")
endif()

# If no build type is specified, default to Debug 
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, defaulting to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# Set the C++ standard
set(HypiC_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD ${HypiC_CXX_STANDARD} CACHE STRING "Default C++ standard")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require C++ standard")

# Set module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Check that the compiler supports the C++ standard
set(HypiC_MIN_GCC_VERSION 8.0)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${HypiC_MIN_GCC_VERSION})
    message(FATAL_ERROR "GCC version must be at least ${HypiC_MIN_GCC_VERSION}")
  endif()
else()
  message(WARNING 
          "You are using an unsupported compiler! "
          "Please use GCC >= ${HypiC_MIN_GCC_VERSION}")
endif()


#add warnings?


# Source files
set(HypiC_SOURCES
"src/Particles_Object.cpp"
"src/Rate_Table.cpp"
)

add_library(HypiC SHARED ${HypiC_SOURCES})

target_include_directories(HypiC
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

#===============================================================================
# Other directories
#===============================================================================

# Enable testing
include(CTest)
enable_testing()
add_subdirectory(test)