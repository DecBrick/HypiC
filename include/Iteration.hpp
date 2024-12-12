#pragma once
#include "Particles_Object.hpp"
#include "math.hpp"

namespace HypiC{

    HypiC::Particles_Object Update_Heavy_Species_Neutrals(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_rates, HypiC::Options_Object Simulation_Parameters);
    HypiC::Particles_Object Update_Heavy_Species_Ions(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_rates, HypiC::Options_Object Simulation_Parameters);

    HypiC::Electrons_Object Update_Electrons(HypiC::Electrons_Object Electrons, HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_Rates, HypiC::Rate_Table_Object Collision_Loss_Rates, HypiC::Options_Object Simulation_Parameters);//, double t);

    //double Freq_Electron_Ion(double EnergyDensity, double ElectronTemp, double IonZ);

    //double Coulomb_Logarithm(double EnergyDensity, double ElectronTemp, double IonZ);

    void Solve_Potential(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters);


    void Compute_Electric_Field(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters, double Discharge_Current);

    double Integrate_Discharge_Current(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters);

    void Compute_Pressure_Gradient(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters);

    void Update_Thermal_Conductivity(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters);

    void Update_Electron_Energy(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters, HypiC::Rate_Table_Object Loss_Rates);
}