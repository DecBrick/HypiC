#include "HypiCpp.hpp"
#include <string>

//testing library
#include "unit_test_framework.hpp"


TEST_CASE(Constructor){
    //make the instance 
    //HypiC::Rate_Table_Object<double> Test_Table = HypiC::Rate_Table_Object<double>();
    HypiC::Rate_Table_Object Test_Table = HypiC::Rate_Table_Object();

    ASSERT(Test_Table._Energies == nullptr);
    ASSERT(Test_Table._Rates == nullptr);


}

TEST_CASE(Xe_Single_Ionization){
    //make the instance 
    //HypiC::Rate_Table_Object<double> Test_Table = HypiC::Rate_Table_Object<double>();
    HypiC::Rate_Table_Object Test_Table = HypiC::Rate_Table_Object();

    //call the read file
    std::string base_dir = "${CMAKE_SOURCE_DIR}";
    std::string file_path = "Reactions/Xe_Ionization_0_to_1.txt";
    Test_Table.Read_Table(base_dir + file_path);

    //test that file was read correctly
    ASSERT(1 == 1);

    //interpolation tests
}

TEST_SUITE(dummy_suite){
    TEST(Constructor);
    TEST(Xe_Single_Ionization);
}



auto 
main() -> int
{
    RUN_SUITE(dummy_suite);
    return 0;
}