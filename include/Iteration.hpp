#pragma once
#include "Particles_Object.hpp"

namespace HypiC{

    void Update_Heavy_Species(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_Rates, HypiC::Options_Object Simulation_Parameters);

    void Update_Electrons(HypiC::Electrons_Object Electrons, HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_Rates, HypiC::Options_Object Simulation_Parameters);//, double t);

    //double Freq_Electron_Ion(double EnergyDensity, double ElectronTemp, double IonZ);

    //double Coulomb_Logarithm(double EnergyDensity, double ElectronTemp, double IonZ);

    void Solve_Potential(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters);

    double Forward_Difference(double f0, double f1, double f2, double x0, double x1, double x2);

    double Central_Difference(double f0, double f1, double f2, double x0, double x1, double x2);

    double Backward_Difference(double f0, double f1, double f2, double x0, double x1, double x2);

    double Linear_Transition(double x, double cutoff, double L, double y1, double y2);

    void Compute_Electric_Field(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters, double Discharge_Current);

    double Integrate_Discharge_Current(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters);

    void Compute_Pressure_Gradient(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters);
}