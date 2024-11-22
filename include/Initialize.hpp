#pragma once 
#include <cstddef>
#include <string>
#include <cmath>
#include "Particles_Object.hpp"
#include "Output.hpp"

namespace HypiC{
    class Options_Object{
        private:

        public:
            //general input options
            size_t nIterations;
            size_t Output_Interval;
            //geometric input options
            double Domain_Length_m;
            size_t nCells;
            double Channel_Length_m;
            //operating conditions
            double Discharge_Voltage_V;
            double Mass_Flow_Rate_kg_s;
            //particle input options
            size_t N_Neutrals;
            size_t N_Ions; 
            double Initial_Neutral_Temperature_K; 
            double Initial_Ion_Temperature_eV; 

            //electron options
            double Initial_Max_Electron_Temperature_eV;
            double Initial_Min_Ion_Density = 2e17;
            double Initial_Max_Ion_Density = 1e18;
            double Initial_Anode_Temperature_eV;
            double Initial_Cathode_Temperature_eV;
    };


    HypiC::Options_Object Read_Input(std::string Filename)
    {
        //construct the object
        HypiC::Options_Object read_options = HypiC::Options_Object();
    };

    //for the initialization of particles, I think we can physically distribute the particles evenly
    //velocities come from maxwellian
    //weights come from density 
    //electron density and electron temperature come from the inital distributions
    //electric field can be set to 0 then updated (update electrons first?). 

    double Initial_Electron_Density(double z, double n_min, double n_max, double Vd, double mdot, double Lch);

    double Initial_Electron_Temperature(double z, double Te_Anode, double Te_Cathode, double Te_Max, double Lch, double z_max);

    double Initial_Ion_Bulk_Velocity(double Te_Anode, double Vd, double z, double Lch, double z_max);

    HypiC::Particles_Object Initialize_Neutrals(HypiC::Options_Object Inputs);

    HypiC::Particles_Object Initialize_Ions(HypiC::Options_Object Inputs);

    HypiC::Time_Sum_Object Zero_Time_Sum()
    {
        
    };
}