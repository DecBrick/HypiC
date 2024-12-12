#pragma once 
#include <cstddef>
#include <string>
#include <cmath>
#include "Particles_Object.hpp"
#include "Electrons_Object.hpp"


namespace HypiC{

    HypiC::Electrons_Object Particles_to_Grid(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons);

    HypiC::Particles_Object Grid_to_Particles_Neutrals(HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons);
    HypiC::Particles_Object Grid_to_Particles_Ions(HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons);
}