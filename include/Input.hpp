#pragma once
#include <string>

namespace HypiC{
    class Options_Object{
        private:

        public:
            //general input options
            size_t nIterations;
            size_t Output_Interval;
            double dt;
            size_t nCells;
            //geometric input options
            double Domain_Length_m;
            double Channel_Length_m;
            double Channel_Area_m2;
            //operating conditions
            double Discharge_Voltage_V;
            double Mass_Flow_Rate_kg_s;
            //particle input options
            size_t N_Neutrals;
            size_t N_Ions;
            double Initial_Neutral_Temperature_K; 
            double Initial_Ion_Temperature_K; 

            //electron options
            double Initial_Max_Electron_Temperature_eV;
            double Initial_Anode_Temperature_eV;
            double Initial_Cathode_Temperature_eV;
            double Initial_Min_Ion_Density;
            double Initial_Max_Ion_Density;
            double Min_Electron_Temperature_eV;


            //construction/destruction methods
            Options_Object();
            ~Options_Object();
            //Read from file 
            void Read_Input(std::string Filename);
    };
}