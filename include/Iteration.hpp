#pragma once
#include "Particles_Object.hpp"

namespace HypiC{

    void Update_Heavy_Species(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_Rates, double dt);

}