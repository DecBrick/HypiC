#include "HypiCpp.hpp"
#include <cmath>

//testing library
#include "unit_test_framework.hpp"

double Tolerance = 1e-12;

TEST_CASE(Construction){
    //initialize
    HypiC::Particles_Object Part = HypiC::Particles_Object();
    //check that the particle counter and properties initialized correctly
    ASSERT(Part._nParticles == 0);
    ASSERT(Part._ChargetoMassRatio == 0);
    ASSERT(Part._IonizationDirection == 0);
    //check that the vectors are empty 
    ASSERT(Part._Positions.empty());
    ASSERT(Part._Velocities.empty());
    ASSERT(Part._Weights.empty());
    ASSERT(Part._ElectricField.empty());
    ASSERT(Part._ElectronDensity.empty());
    ASSERT(Part._ElectronTemperature.empty());

}

TEST_CASE(Add_Remove_Particles){
    //initialize
    HypiC::Particles_Object Part = HypiC::Particles_Object();

    //check that we can add particles 
    Part.Add_Particle(1.0, 1000, 1, 10000, 1e18, 10);

    //pull properties
    //also is testing accessor methods
    ASSERT(Part._nParticles == 1);
    ASSERT_NEAR(Part.get_Position(0), 1.0, Tolerance);
    ASSERT_NEAR(Part.get_Velocity(0), 1000.0, Tolerance);
    ASSERT_NEAR(Part.get_Weight(0), 1.0, Tolerance);
    ASSERT_NEAR(Part.get_ElectricField(0), 10000.0, Tolerance);
    ASSERT_NEAR(Part.get_ElectronDensity(0), 1e18, Tolerance);
    ASSERT_NEAR(Part.get_ElectronTemperature(0), 10.0, Tolerance);

    //check that we can remove particles
    Part.Remove_Particle(0);
    //check that the values revert to what we expect
    ASSERT(Part._nParticles == 0);
    //check that the vectors are empty 
    ASSERT(Part._Positions.empty());
    ASSERT(Part._Velocities.empty());
    ASSERT(Part._Weights.empty());
    ASSERT(Part._ElectricField.empty());
    ASSERT(Part._ElectronDensity.empty());
    ASSERT(Part._ElectronTemperature.empty());
}

TEST_CASE(Update_Particles){
    //initialize and add particle
    HypiC::Particles_Object Ion_Test = HypiC::Particles_Object();
    HypiC::Particles_Object Neutral_Test = HypiC::Particles_Object();
    Neutral_Test.Add_Particle(1.0, 1000, 1, 10000, 1e18, 10);
    Ion_Test.Add_Particle(1.0, 1000, 1, 10000, 1e18, 10);

    //set ionization directions
    Neutral_Test._IonizationDirection = -1;
    Ion_Test._IonizationDirection = 1;

    //initialize ionization rates
    HypiC::Rate_Table_Object Test_Table = HypiC::Rate_Table_Object();
    //this assumes that the build directory is in the same folder as the HypiC directory, should enforce this.
    std::string file_path = "../../HypiC/Reactions/Xe_Ionization_0_to_1.txt";
    Test_Table.Read_Table(file_path);

    //update the neutral particle 
    Neutral_Test.Update_Particle(0, 1e-9, Test_Table);

    //slight change to position, velocities unchanged, weight updated by neutralization 
    ASSERT_NEAR(Neutral_Test.get_Position(0), 1.0 + 1e-6, Tolerance);
    ASSERT_NEAR(Neutral_Test.get_Velocity(0), 1000.0, Tolerance);
    ASSERT_NEAR(Neutral_Test.get_Weight(0), exp(-1 * 1e18 * 1e-9 * 1.655e-14), Tolerance);

    //set Xe charge to mass ratio 
    Ion_Test._ChargetoMassRatio = 1.6e-19 / (131.29 * 1.67e-27);
    //Update the ion 
    Ion_Test.Update_Particle(0, 1e-9, Test_Table);

    //velocity slightly increased 
    ASSERT_NEAR(Ion_Test.get_Velocity(0), 1000 + 1.6e-24 / (131.29 * 1.67e-27), Tolerance);
    ASSERT_NEAR(Ion_Test.get_Position(0), 1 + 1e-9 *(1000 + 1.6e-24 / (131.29 * 1.67e-27)), Tolerance);
    ASSERT_NEAR(Ion_Test.get_Weight(0),  exp(1 * 1e18 * 1e-9 * 1.655e-14), Tolerance);

    //check set methods
    Ion_Test.set_Position(0, 2.0);
    Ion_Test.set_Velocity(0, 500.0);

    ASSERT_NEAR(Ion_Test.get_Position(0), 2.0, Tolerance);
    ASSERT_NEAR(Ion_Test.get_Velocity(0), 500.0, Tolerance);
}

TEST_SUITE(Particles_Suite){
    TEST(Construction);
    TEST(Add_Remove_Particles);
    TEST(Update_Particles);
}



auto 
main() -> int
{
    RUN_SUITE(Particles_Suite);
    return 0;
}