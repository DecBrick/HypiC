#include <string>

#include "HypiCpp.hpp"
#include <iostream>


int main(){
    std::cout << "-------------------------------------\n";
    std::cout << "HypiC++ version 0.0.0\n";
    std::cout << "-------------------------------------\n";
    //read input file, should be caled Input in the same folder as the executable 
    HypiC::Options_Object Input_Options = HypiC::Options_Object();
    Input_Options.Read_Input("Input.txt");

    std::cout << "Input File Read\n";

    //read the ionization and loss rates. 
    std::cout << "Initializating Collisional Rates\n";
    HypiC::Rate_Table_Object Ionization_Rates = HypiC::Rate_Table_Object();
    HypiC::Rate_Table_Object Loss_Rates = HypiC::Rate_Table_Object();
    //this assumes that the build directory is in the same folder as the HypiC directory, should enforce this.
    std::string file_path1 = "../HypiC/Reactions/Xe_Ionization_0_to_1.txt";
    Ionization_Rates.Read_Table(file_path1);
    std::string file_path2 = "../HypiC/Reactions/Xe_Loss.txt";
    Loss_Rates.Read_Table(file_path2);

    std::cout << "Reaction Rates Read\n";

    //initialize objects
    HypiC::Particles_Object Neutrals = HypiC::Initialize_Neutrals(Input_Options);
    HypiC::Particles_Object Ions = HypiC::Initialize_Ions(Input_Options);
    HypiC::Electrons_Object Electrons = HypiC::Initialize_Electrons(Input_Options);

    std::cout << "Objects\n";
    //initial interpolations 
    //interpolate particle densities to electrons
    HypiC::Particles_to_Grid(Neutrals, Ions, Electrons);

    //solve electric field for electrons, add ji to the electrons object 
    HypiC::Update_Mobility(Electrons, Input_Options, Ionization_Rates);
    HypiC::Compute_Pressure_Gradient(Electrons, Input_Options);
    double Discharge_Current = HypiC::Integrate_Discharge_Current(Electrons, Input_Options);
    HypiC::Compute_Electric_Field(Electrons, Input_Options, Discharge_Current); //or whatever the function is
    //interpolate field back to particles
    HypiC::Grid_to_Particles(Neutrals, Ions, Electrons);

    //take back half step for ions (neutrals are unaffected by the field) 
    Ions.Velocity_Backstep(Input_Options.dt);
    std::cout << "Plasma Objects Initialized\n";

    //initialize outputs 
    std::cout << "Initializating Output\n";
    HypiC::Time_Sum_Object Results = HypiC::Time_Sum_Object();
    Results.Initialize_Time_Sum(Input_Options.nCells, Electrons);
    
    

    std::cout << "Initialization Complete\n";

    //main loop
    for(size_t i=0; i < Input_Options.nIterations; ++i){
        //update heavy species
        HypiC::Update_Heavy_Species(Neutrals, Ions, Ionization_Rates, Input_Options);

        //interpolate
        HypiC::Particles_to_Grid(Neutrals, Ions, Electrons);

        //update electrons
        HypiC::Update_Electrons(Electrons, Neutrals, Ions, Ionization_Rates, Loss_Rates, Input_Options);

        //interpolate
        HypiC::Grid_to_Particles(Neutrals, Ions, Electrons);

        //update time sum
        Results.Time_Sum(Electrons, Input_Options);

        //periodic output
        if (i % Input_Options.Output_Interval == 0){
            Results.Write_Output("Output.csv", Input_Options.nCells);
        }
        
    }

    //final output
    std::cout << "Simulation Complete\n";
    Results.Write_Output("Output.csv", Input_Options.nCells);
}