include(${PROJECT_SOURCE_DIR}/cmake/HypiC_add_test.cmake)

#add test files as fit
HypiC_add_test(./Particle_Tests.cpp)
HypiC_add_test(./Rate_Table_Tests.cpp)
HypiC_add_test(./Initialization_Tests.cpp)
HypiC_add_test(./Output_Tests.cpp)
HypiC_add_test(./Math_Tests.cpp)