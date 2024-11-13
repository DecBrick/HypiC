#include <HypiCpp.hpp>


//testing library
#include "unit_test_framework.hpp"


TEST_CASE(dummy){
    ASSERT( 1 == 1);
    HypiC::Particles_Object<double> Part = HypiC::Particles_Object<double>();
}

TEST_SUITE(dummy_suite){
    TEST(dummy);
}



auto 
main() -> int
{
    RUN_SUITE(dummy_suite);
    return 0;
}