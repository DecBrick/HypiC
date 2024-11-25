#pragma once 
#include <cstddef>
#include <string>
#include <cmath>
#include "Particles_Object.hpp"
#include "Electrons_Object.hpp"
#include "Output.hpp"

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


            //construction/destruction methods
            Options_Object();
            ~Options_Object();
            //Read from file 
            void Read_Input(std::string Filename);
    };
    

    double Initial_Electron_Density(double z, double n_min, double n_max, double Vd, double mdot, double Lch);

    double Initial_Electron_Temperature(double z, double Te_Anode, double Te_Cathode, double Te_Max, double Lch, double z_max);

    double Initial_Ion_Bulk_Velocity(double Te_Anode, double Vd, double z, double Lch, double z_max);

    HypiC::Particles_Object Initialize_Neutrals(HypiC::Options_Object Inputs);

    HypiC::Particles_Object Initialize_Ions(HypiC::Options_Object Inputs);

    HypiC::Electrons_Object Initialize_Electrons(HypiC::Options_Object Inputs);

    //HypiC::Time_Sum_Object Zero_Time_Sum();
    

    double Maxwellian_Sampler(double mu, double sigma);
}