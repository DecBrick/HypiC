#pragma once 
#include <cstddef>
#include <string>
#include <cmath>
#include "Input.hpp"
#include "Particles_Object.hpp"
#include "Electrons_Object.hpp"


namespace HypiC{
    
    double Initial_Magnetic_Field(double B_max, double Lch, double z);

    double Initial_Electron_Density(double z, double n_min, double n_max, double Vd, double mdot, double Lch);

    double Initial_Electron_Temperature(double z, double Te_Anode, double Te_Cathode, double Te_Max, double Lch, double z_max);

    double Initial_Ion_Bulk_Velocity(double Te_Anode, double Vd, double z, double Lch, double z_max);

    HypiC::Particles_Object Initialize_Neutrals(HypiC::Options_Object Inputs);

    HypiC::Particles_Object Initialize_Ions(HypiC::Options_Object Inputs);

    HypiC::Electrons_Object Initialize_Electrons(HypiC::Options_Object Inputs);
    

    double Maxwellian_Sampler(double mu, double sigma);
}