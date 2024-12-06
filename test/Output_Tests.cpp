#include "HypiCpp.hpp"
#include <cmath>

//testing library
#include "unit_test_framework.hpp"


double Tolerance = 1e-12;

TEST_CASE(Write_File){
    double z;
    //make the object
    HypiC::Time_Sum_Object T_sum = HypiC::Time_Sum_Object();

    HypiC::Electrons_Object e = HypiC::Electrons_Object();

    double dz = 0.05 / 200;
    for(size_t c=0; c<200; ++c){
        z = 0.5 * dz + dz * c;
        e.Cell_Center.push_back(z);
        //std::cout << e.Cell_Center[c] << "\n";
    }

    //call the set values
    T_sum.Initialize_Time_Sum(200, e);


    //write the file to the current directory
    T_sum.Write_Output("Test_Output.csv", 200);

    //tell the user to look for the file. 
    std::cout << "Please check the testing directory for a file called Test_Output.csv";
}


TEST_SUITE(Output_Suite){
    TEST(Write_File);
}

auto 
main() -> int
{
    RUN_SUITE(Output_Suite);
    return 0;
}