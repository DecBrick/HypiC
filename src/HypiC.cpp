#include <string>

#include "HypiCpp.hpp"
#include <iostream>
#include <chrono>


int main(){
    auto start = std::chrono::high_resolution_clock::now();
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
    Electrons = HypiC::Particles_to_Grid(Neutrals, Ions, Electrons);

    //solve electric field for electrons, add ji to the electrons object
    Electrons.Update_Mobility(Input_Options, Ionization_Rates);
    Electrons.Update_Pressure_Gradient(Input_Options);
    Electrons.Compute_Discharge_Current(Input_Options);
    Electrons.Compute_Electric_Field(Input_Options); 

    std::cout << "Electrons \n";
    //interpolate field back to particles
    Neutrals = HypiC::Grid_to_Particles_Neutrals(Neutrals, Electrons);
    Ions = HypiC::Grid_to_Particles_Ions(Ions, Electrons);

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
        //std::cout << "Neutral Updated\n";
        //Ions = HypiC::Update_Heavy_Species_Ions(Neutrals, Ions, Ionization_Rates, Input_Options);
        //std::cout << "Ions Updated\n";
        //interpolate
        std::cout << "Updated\n";
        Electrons = HypiC::Particles_to_Grid(Neutrals, Ions, Electrons);

        //update electrons
        Electrons = HypiC::Update_Electrons(Electrons, Neutrals, Ions, Ionization_Rates, Loss_Rates, Input_Options);

        //interpolate
        Neutrals = HypiC::Grid_to_Particles_Neutrals(Neutrals, Electrons);
        Ions = HypiC::Grid_to_Particles_Ions(Ions, Electrons);

        //update time sum
        Results.Time_Sum(Electrons, Input_Options);

        //periodic output
        if (i % Input_Options.Output_Interval == 0){
            Results.Write_Output("Output.csv", Input_Options.nCells);
        }
        
    }

    //final output
    std::cout << "Simulation Complete\n";
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);

    std::cout << "Run time in s was:" << duration.count() << "\n";

    Results.Write_Output("Output.csv", Input_Options.nCells);
}