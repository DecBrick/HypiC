#include "HypiCpp.hpp"
#include <string>
#include <vector>

//testing library
#include "unit_test_framework.hpp"

double Tolerance = 1e-12;

TEST_CASE(Thomas_Algorithm){
    std::vector<double> ld;
    std::vector<double> d;
    std::vector<double> ud;
    std::vector<double> B; 
    std::vector<double> x; 

    //add lower diagonal
    ld.push_back(5.5);
    ld.push_back(11);
    ld.push_back(M_PI);
    ld.push_back(0.0);

    //add diagonal 
    d.push_back(1.0);
    d.push_back(3.0);
    d.push_back(9.1);
    d.push_back(2.5);
    d.push_back(M_PI * M_PI);

    //add upper diagonal
    ud.push_back(2.0);
    ud.push_back(4.5);
    ud.push_back(20.0);
    ud.push_back(2.71);


    //add B, use manufactured solution where x is a vector of 1's
    B.push_back(3.0);
    B.push_back(13.0);
    B.push_back(40.1);
    B.push_back(M_PI + 5.21);
    B.push_back(M_PI * M_PI);

    x = HypiC::Thomas_Algorithm(ld, d, ud, B);

    ASSERT_NEAR(x[0], 1, Tolerance);
    ASSERT_NEAR(x[1], 1, Tolerance);
    ASSERT_NEAR(x[2], 1, Tolerance);
    ASSERT_NEAR(x[3], 1, Tolerance);
    ASSERT_NEAR(x[4], 1, Tolerance);
}


TEST_SUITE(Math){
    TEST(Thomas_Algorithm);
}



auto 
main() -> int
{
    RUN_SUITE(Math);
    return 0;
}