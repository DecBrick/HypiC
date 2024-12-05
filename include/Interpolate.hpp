#pragma once 
#include <cstddef>
#include <string>
#include <cmath>
#include "Particles_Object.hpp"
#include "Electrons_Object.hpp"


namespace HypiC{

    void Particles_to_Grid(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons);

    void Grid_to_Particles(HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons);
}