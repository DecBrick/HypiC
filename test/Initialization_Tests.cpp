#include "HypiCpp.hpp"
#include <cmath>
#include <string>
#include <algorithm>
#include <iomanip>

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
    ASSERT_NEAR(Inputs.Domain_Length_m, 0.05, Tolerance);
    ASSERT_NEAR(Inputs.Channel_Length_m, 0.025 ,Tolerance);
    ASSERT_NEAR(Inputs.Channel_Area_m2, M_PI * (pow(0.05, 2)-pow(0.035, 2)) ,Tolerance);
    ASSERT_NEAR(Inputs.Discharge_Voltage_V, 300, Tolerance);
    ASSERT_NEAR(Inputs.Mass_Flow_Rate_kg_s, 5e-6, Tolerance);
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

TEST_CASE(Initialize_Neutrals){
    //read the example input options
    //this assumes that the build directory is in the same folder as the HypiC directory, should enforce this.
    std::string file_path = "../../HypiC/Example_Input.txt";
    HypiC::Options_Object Inputs = HypiC::Options_Object();
    Inputs.Read_Input(file_path);

    //call the initialize function
    HypiC::Particles_Object Neutrals = HypiC::Initialize_Neutrals(Inputs);

    //check that the number of particles is correct
    ASSERT(Neutrals._nParticles == Inputs.N_Neutrals);

    //check for the charge to mass ratio and ionization direction
    ASSERT_NEAR(Neutrals._ChargetoMassRatio, 0, Tolerance);

    //check that positions are within bounds

    ASSERT(*std::max_element(Neutrals._Positions.begin(), Neutrals._Positions.end()) <= Inputs.Domain_Length_m);
    ASSERT(*std::min_element(Neutrals._Positions.begin(), Neutrals._Positions.end()) >= 0.0);
    //check that velocity is physical
    //I don't think anything should be going faster than 1000 km/s
    ASSERT(*std::max_element(Neutrals._Velocities.begin(), Neutrals._Velocities.end()) <= 1000000);
    ASSERT(*std::min_element(Neutrals._Velocities.begin(), Neutrals._Velocities.end()) >= -1000000);

}

TEST_CASE(Initialize_Ions){ 
    //read the example input options
    //this assumes that the build directory is in the same folder as the HypiC directory, should enforce this.
    std::string file_path = "../../HypiC/Example_Input.txt";
    HypiC::Options_Object Inputs = HypiC::Options_Object();
    Inputs.Read_Input(file_path);

    //call the initialize function
    HypiC::Particles_Object Ions = HypiC::Initialize_Ions(Inputs);

    //check that the number of particles is correct
    ASSERT(Ions._nParticles == Inputs.N_Ions);

    //check for the charge to mass ratio and ionization direction
    double mass_charge = 1.602176634e-19 / (131.29 * 1.66053907e-27);
    ASSERT_NEAR(Ions._ChargetoMassRatio, mass_charge, Tolerance);

    //check that positions are within bounds

    ASSERT(*std::max_element(Ions._Positions.begin(), Ions._Positions.end()) <= Inputs.Domain_Length_m);
    ASSERT(*std::min_element(Ions._Positions.begin(), Ions._Positions.end()) >= 0.0);
    //check that velocity is physical
    //I don't think anything should be going faster than 1000 km/s
    ASSERT(*std::max_element(Ions._Velocities.begin(), Ions._Velocities.end()) <= 1000000);
    ASSERT(*std::min_element(Ions._Velocities.begin(), Ions._Velocities.end()) >= -1000000);

}

TEST_SUITE(Initalization){
    TEST(Read_Input);
    TEST(Initialize_Neutrals);
    TEST(Initialize_Ions)
}


auto 
main() -> int
{
    RUN_SUITE(Initalization);
    return 0;
}