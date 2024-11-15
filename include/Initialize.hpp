#pragma once 
#include <cstddef>
#include <string>
#include "Particles_Object.hpp"
#include "Output.hpp"

namespace HypiC{
    class Options_Object{
        private:

        public:
            size_t nIterations;
            size_t Output_Interval;
    };


    HypiC::Options_Object Read_Input(std::string Filename)
    {
        
    };

    HypiC::Particles_Object Initialize_Neutrals(HypiC::Options_Object Inputs)
    {
        
    };

    HypiC::Particles_Object Initialize_Ions(HypiC::Options_Object Inputs)
    {
        
    };

    HypiC::Time_Sum_Object Zero_Time_Sum()
    {
        
    };
}