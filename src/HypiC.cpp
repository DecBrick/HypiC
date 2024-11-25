#include <string>

#include "HypiCpp.hpp"


void main(){

    //read input file, should be caled Input in the same folder as the executable 
    HypiC::Options_Object Input_Options = HypiC::Options_Object();
    Input_Options.Read_Input("Input.txt");

    //initialize objects
    HypiC::Particles_Object Neutrals = HypiC::Initialize_Neutrals(Input_Options);
    HypiC::Particles_Object Ions = HypiC::Initialize_Ions(Input_Options);
    //Electrons = HypiC::Initialize_Electrons(Input_Options);

    //initialize outputs 
    HypiC::Time_Sum_Object Results = HypiC::Time_Sum_Object();
    Results.Initialize_Time_Sum(Input_Options.nCells);
    
    //read the ionization rates. 
    HypiC::Rate_Table_Object Ionization_Rates = HypiC::Rate_Table_Object();
    //this assumes that the build directory is in the same folder as the HypiC directory, should enforce this.
    std::string file_path = "../../HypiC/Reactions/Xe_Ionization_0_to_1.txt";
    Ionization_Rates.Read_Table(file_path);

    //main loop
    for(size_t i=0; i < Input_Options.nIterations; ++i){
        //update heavy species
        HypiC::Update_Heavy_Species(Neutrals, Ions, Ionization_Rates, Input_Options);

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
            Results.Write_Output("Output.csv", Input_Options.nCells);
        }
        
    }

    //final output
    Results.Write_Output("Output.csv", Input_Options.nCells);
}