#include <string>

#include "HypiCpp.hpp"


void main(){

    //read input file, should be caled Input in the same folder as the executable 
    HypiC::Options_Object Input_Options = HypiC::Options_Object();
    Input_Options.Read_Input("Input.txt");

    //initialize objects
    HypiC::Particles_Object Neutrals = HypiC::Initialize_Neutrals(Input_Options);
    HypiC::Particles_Object Ions = HypiC::Initialize_Ions(Input_Options);
    HypiC::Electrons_Object Electrons = HypiC::Initialize_Electrons(Input_Options);

    //initial interpolations 
    //interpolate particle densities to electrons
    HypiC::Particles_to_Grid(Neutrals, Ions, Electrons);

    //solve electric field for electrons, add ji to the electrons object 
    double Discharge_Current = Integrate_Discharge_Current(Electrons, Input_Options);
    HypiC::Compute_Electric_Field(Electrons, Input_Options, Discharge_Current); //or whatever the function is

    //interpolate field back to particles
    HypiC::Grid_to_Particles(Ions, Electrons);

    //take back half step for ions (neutrals are unaffected by the field) 
    Ions.Velocity_Backstep(Input_Options.dt);

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
        HypiC::Particles_to_Grid(Neutrals, Ions, Electrons);

        //update electrons
        HypiC::Update_Electrons(Electrons, Neutrals, Ions, Ionization_Rates, Input_Options);

        //interpolate
        HypiC::Grid_to_Particles(Ions, Electrons);

        //update time sum
        Results.Time_Sum(Electrons, Input_Options); 

        //periodic output
        if (i % Input_Options.Output_Interval == 0){
            Results.Write_Output("Output.csv", Input_Options.nCells);
        }
        
    }

    //final output
    Results.Write_Output("Output.csv", Input_Options.nCells);
}