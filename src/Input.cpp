#include <fstream>
#include "HypiCpp.hpp"
#include <iostream>

namespace HypiC{
    //class functions
    Options_Object::Options_Object()
    {
        
    }
    
    // Destructor
    //template <class fp_type>
    //Particles_Object<fp_type>::~Particles_Object()
    Options_Object::~Options_Object()
    {
        
    }
    void Options_Object::Read_Input(std::string Filename){
        std::string line;
        size_t pos;
    
        //open the file and read line by line
        //open the file 
        std::ifstream f(Filename);
        
        //grab line, allow for comments in the input file, each line should end with //
        std::getline(f, line);
        pos = line.find_first_of("//");
        //line 1 is nIterations
        this->nIterations = std::stoul(line.substr(0, pos));
        //rinse and repeat, next line is Output_Interval
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Output_Interval = std::stoul(line.substr(0, pos));
        //dt
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->dt = std::stod(line.substr(0, pos));
        //nCells
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->nCells = std::stoul(line.substr(0, pos));
        //Domain Length
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Domain_Length_m = std::stod(line.substr(0, pos));
        //Channel Length
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Channel_Length_m = std::stod(line.substr(0, pos));
        //Channel Area
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Channel_Area_m2 = std::stod(line.substr(0, pos));
        //Discharge Voltage
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Discharge_Voltage_V = std::stod(line.substr(0, pos));
        //Mass Flow Rate
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Mass_Flow_Rate_kg_s = std::stod(line.substr(0, pos));
        //Number of neutral particles
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->N_Neutrals = std::stoul(line.substr(0, pos));
        //Number of ion particles
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->N_Ions = std::stoul(line.substr(0, pos));
        //Neutral Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Neutral_Temperature_K = std::stod(line.substr(0, pos));
        //Ion Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Ion_Temperature_K = std::stod(line.substr(0, pos));
        //Electron Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Max_Electron_Temperature_eV = std::stod(line.substr(0, pos));
        //Anode Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Anode_Temperature_eV = std::stod(line.substr(0, pos));
        //Cathode Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Cathode_Temperature_eV = std::stod(line.substr(0, pos));
        //Initial_Min_Ion_Denisty
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Min_Ion_Density = std::stod(line.substr(0, pos));
        //Max Ion Density
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Max_Ion_Density = std::stod(line.substr(0, pos));
        //Min Electron Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Min_Electron_Temperature_eV = std::stod(line.substr(0, pos));
    };
}