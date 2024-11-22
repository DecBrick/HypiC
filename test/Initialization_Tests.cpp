#include "HypiCpp.hpp"
#include <cmath>
#include <string>


//testing library
#include "unit_test_framework.hpp"

double Tolerance = 1e-12;

TEST_CASE(Read_Input){
    //create an instance
    HypiC::Options_Object Inputs = HypiC::Options_Object();

    //call the read file
    //this assumes that the build directory is in the same folder as the HypiC directory, should enforce this.
    std::string file_path = "../../HypiC/Example_Input.txt";
    Inputs.Read_Input(file_path);

    //now check that the options are valid
    ASSERT(Inputs.nIterations == 10000);
    ASSERT(Inputs.Output_Interval == 1000);
    ASSERT_NEAR(Inputs.dt, 15e-9 ,Tolerance);
    ASSERT(Inputs.nCells == 200);
    ASSERT_NEAR(Inputs.Domain_Length_m, 0.1, Tolerance);
    ASSERT_NEAR(Inputs.Channel_Length_m, 0.025 ,Tolerance);
    ASSERT_NEAR(Inputs.Discharge_Voltage_V, 300, Tolerance);
    ASSERT_NEAR(Inputs.Mass_Flow_Rate_kg_s, 5.16e-6, Tolerance);
    ASSERT(Inputs.N_Neutrals == 10000);
    ASSERT(Inputs.N_Ions == 10000);
    ASSERT_NEAR(Inputs.Initial_Neutral_Temperature_K, 773, Tolerance);
    ASSERT_NEAR(Inputs.Initial_Ion_Temperature_K, 1160.4, Tolerance);
    ASSERT_NEAR(Inputs.Initial_Max_Electron_Temperature_eV, 30, Tolerance);
    ASSERT_NEAR(Inputs.Initial_Anode_Temperature_eV, 3, Tolerance);
    ASSERT_NEAR(Inputs.Initial_Cathode_Temperature_eV, 5, Tolerance);
    ASSERT_NEAR(Inputs.Initial_Min_Ion_Density, 2e17, Tolerance);
    ASSERT_NEAR(Inputs.Initial_Max_Ion_Density, 1e18, Tolerance);
}

TEST_CASE(Aux_Functions){

}

TEST_CASE(Initialize_Neutrals){

}

TEST_CASE(Initialize_Ions){ 

}

TEST_SUITE(Initalization){
    TEST(Read_Input);
    TEST(Aux_Functions);
    TEST(Initialize_Neutrals);
    TEST(Initialize_Ions)
}


auto 
main() -> int
{
    RUN_SUITE(Initalization);
    return 0;
}