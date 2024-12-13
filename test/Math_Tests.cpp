#include "HypiCpp.hpp"
#include <string>

//testing library
#include "unit_test_framework.hpp"


TEST_CASE(Constructor){
    //make the instance 
    //HypiC::Rate_Table_Object<double> Test_Table = HypiC::Rate_Table_Object<double>();
    HypiC::Rate_Table_Object Test_Table = HypiC::Rate_Table_Object();
}
//reading the ionization cross sections
TEST_CASE(Xe_Single_Ionization_Read){
    //make the instance 
    //HypiC::Rate_Table_Object<double> Test_Table = HypiC::Rate_Table_Object<double>();
    HypiC::Rate_Table_Object Test_Table = HypiC::Rate_Table_Object();

    //call the read file
    //this assumes that the build directory is in the same folder as the HypiC directory, should enforce this.
    std::string file_path = "../Reactions/Xe_Ionization_0_to_1.txt";
    Test_Table.Read_Table(file_path);
    
    //test that file was read correctly
    //check first and last energy
    ASSERT_NEAR(Test_Table._Energies[0], 0, 1e-12);
    ASSERT_NEAR(Test_Table._Energies[150], 150, 1e-12);
    //check first and last rate
    ASSERT_NEAR(Test_Table._Rates[0], 0, 1e-12);
    ASSERT_NEAR(Test_Table._Rates[150], 2.983e-13, 1e-14);

}
//Interpolation on cross sections
TEST_CASE(Interpolation){
    //make the instance 
    //HypiC::Rate_Table_Object<double> Test_Table = HypiC::Rate_Table_Object<double>();
    HypiC::Rate_Table_Object Test_Table = HypiC::Rate_Table_Object();

    //call the read file
    //this assumes that the build directory is in the same folder as the HypiC directory, should enforce this.
    std::string file_path = "../Reactions/Xe_Ionization_0_to_1.txt";
    Test_Table.Read_Table(file_path);

    //interpolation tests
    //check bounds
    ASSERT_NEAR(Test_Table.interpolate(-1), 0, 1e-14);
    ASSERT_NEAR(Test_Table.interpolate(200), 2.983e-13, 1e-14);

    //check general interpolation
    ASSERT_NEAR(Test_Table.interpolate(116.5), 2.730e-13, 1e-14);

    //check interpolation at exact table location
    ASSERT_NEAR(Test_Table.interpolate(100), 2.539e-13, 1e-14);
}

TEST_SUITE(Rate_Suite){
    TEST(Constructor);
    TEST(Xe_Single_Ionization_Read);
    TEST(Interpolation);
}



auto 
main() -> int
{
    RUN_SUITE(Rate_Suite);
    return 0;
}