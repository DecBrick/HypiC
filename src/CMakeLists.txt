# Source files
set(SOURCES
"Electrons_Object.cpp"
"Initialize.cpp"
"Interpolate.cpp"
"Iteration.cpp"
"Output.cpp"
"Particles_Object.cpp"
"Rate_Table.cpp"
"math.cpp"
"Input.cpp"
)

add_library(HypiC_Sources SHARED ${SOURCES})

target_include_directories(HypiC_Sources
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)