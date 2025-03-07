#pragma once
#include "Particles_Object.hpp"
#include "math.hpp"

namespace HypiC{

    HypiC::Particles_Object Update_Heavy_Species_Neutrals(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Grid, HypiC::Options_Object Simulation_Parameters);
    HypiC::Particles_Object Update_Heavy_Species_Ions(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Grid, HypiC::Options_Object Simulation_Parameters);

    HypiC::Electrons_Object Update_Electrons(HypiC::Electrons_Object Electrons, HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_Rates, HypiC::Rate_Table_Object Collision_Loss_Rates, HypiC::Options_Object Simulation_Parameters);//, double t);

    //double Freq_Electron_Ion(double EnergyDensity, double ElectronTemp, double IonZ);

    //double Coulomb_Logarithm(double EnergyDensity, double ElectronTemp, double IonZ);

}