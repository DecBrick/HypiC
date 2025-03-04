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
    //check that the vectors are empty 
    ASSERT(Part._Positions.empty());
    ASSERT(Part._Velocities.empty());
    ASSERT(Part._Weights.empty());
    ASSERT(Part._CellIndex.empty());
    ASSERT(Part._ElectricField.empty());

}

TEST_CASE(Add_Remove_Particles){
    //initialize
    HypiC::Particles_Object Part = HypiC::Particles_Object();

    //check that we can add particles 
    Part.Add_Particle(1.0, 1000, 1, 0, 10000);

    //pull properties
    //also is testing accessor methods
    ASSERT(Part._nParticles == 1);
    ASSERT_NEAR(Part.get_Position(0), 1.0, Tolerance);
    ASSERT_NEAR(Part.get_Velocity(0), 1000.0, Tolerance);
    ASSERT_NEAR(Part.get_Weight(0), 1.0, Tolerance);
    ASSERT(Part.get_Cell(0) == 0);
    ASSERT_NEAR(Part.get_ElectricField(0), 10000.0, Tolerance);

    //check that we can remove particles
    Part.Remove_Particle(0);
    //check that the values revert to what we expect
    ASSERT(Part._nParticles == 0);
    //check that the vectors are empty 
    ASSERT(Part._Positions.empty());
    ASSERT(Part._Velocities.empty());
    ASSERT(Part._Weights.empty());
    ASSERT(Part._CellIndex.empty());
    ASSERT(Part._ElectricField.empty());
}

TEST_CASE(Update_Particles){
    //initialize and add particle
    HypiC::Particles_Object Ion_Test = HypiC::Particles_Object();
    HypiC::Particles_Object Neutral_Test = HypiC::Particles_Object();
    Neutral_Test.Add_Particle(0.6, 1000, 1, 1, 10000);
    Ion_Test.Add_Particle(0.6, 1000, 1, 1, 10000);

    //initialize a basic grid [0.25, 0.75, 1.25]
    HypiC::Electrons_Object Grid = HypiC::Electrons_Object();
    Grid.Add_Electron(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25);
    Grid.Add_Electron(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75);
    Grid.Add_Electron(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.25);
    Grid.Grid_Step = 0.5;

    //update the neutral particle 
    Neutral_Test.Push_Particle(0, 1e-9, Grid);

    //slight change to position, velocities unchanged, cell unchanged 
    ASSERT_NEAR(Neutral_Test.get_Position(0), 0.6 + 1e-6, Tolerance);
    ASSERT_NEAR(Neutral_Test.get_Velocity(0), 1000.0, Tolerance);
    ASSERT(Neutral_Test.get_Cell(0) == 1);


    //add three more neutrals to test the cell index update
    Neutral_Test.Add_Particle(0.49, 1e3, 1000, 0, 0.0);
    Neutral_Test.Add_Particle(0.51, -1e3, 1000, 1, 0.0);
    Neutral_Test.Add_Particle(0.49, 5e4, 1000, 0, 0.0);

    //update these particles
    Neutral_Test.Push_Particle(1, 2e-5, Grid);
    Neutral_Test.Push_Particle(2, 2e-5, Grid);
    Neutral_Test.Push_Particle(3, 2e-5, Grid);

    //check that the cells have been updated correctly
    ASSERT(Neutral_Test.get_Cell(1) == 1);
    ASSERT(Neutral_Test.get_Cell(2) == 0);
    ASSERT(Neutral_Test.get_Cell(3) == 2);

    //set Xe charge to mass ratio 
    Ion_Test._ChargetoMassRatio = 1.6e-19 / (131.29 * 1.67e-27);
    //Update the ion 
    Ion_Test.Push_Particle(0, 1e-9, Grid);

    //velocity slightly increased 
    ASSERT_NEAR(Ion_Test.get_Velocity(0), 1000 + 1.6e-24 / (131.29 * 1.67e-27), Tolerance);
    ASSERT_NEAR(Ion_Test.get_Position(0), 0.6 + 1e-9 *(1000 + 1.6e-24 / (131.29 * 1.67e-27)), Tolerance);
    ASSERT(Ion_Test.get_Cell(0) == 1);

    //check set methods
    Ion_Test.set_Position(0, 2.0);
    Ion_Test.set_Velocity(0, 500.0);

    ASSERT_NEAR(Ion_Test.get_Position(0), 2.0, Tolerance);
    ASSERT_NEAR(Ion_Test.get_Velocity(0), 500.0, Tolerance);

    //check initialization
    Ion_Test.Velocity_Backstep(1e-9);

    ASSERT_NEAR(Ion_Test.get_Velocity(0), 500 - (1.6e-19 / (131.29 * 1.67e-27)* 10000 * 5e-10),Tolerance);
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