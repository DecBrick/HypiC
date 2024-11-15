#include <string>

#include "HypiCpp.hpp"


void main(){

    //read input file
    HypiC::Options_Object Input_Options = HypiC::Read_Input("Input.txt");

    //initialize
    HypiC::Particles_Object Neutrals = HypiC::Initialize_Neutrals(Input_Options);
    HypiC::Particles_Object Ions = HypiC::Initialize_Ions(Input_Options);
    //Electrons = HypiC::Initialize_Electrons(Input_Options);
    HypiC::Time_Sum_Object Results = HypiC::Zero_Time_Sum();


    //main loop
    for(size_t i=0; i < Input_Options.nIterations; ++i){
        //update heavy species
        HypiC::Update_Heavy_Species(Neutrals, Ions);

        //interpolate
        //HypiC::Particles_to_Grid(Neutrals, Ions, Electrons);

        //update electrons
        //HypiC::Update_Electrons(Electrons);

        //interpolate
        //HypiC::Grid_to_Particles(Neutrals, Ions, Electrons);

        //update time sum
        //Results = HypiC::Time_Sum(Neutrals, Ions, Electrons); 

        //periodic output
        if (i % Input_Options.Output_Interval == 0){
            HypiC::Output(Results);
        }
        
    }

    //final output
    HypiC::Output(Results);
}