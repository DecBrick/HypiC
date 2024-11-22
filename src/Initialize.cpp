#include <cmath>
#include "HypiCpp.hpp"

namespace HypiC{
    double Initial_Electron_Density(double z, double n_min, double n_max, double Vd, double mdot, double Lch){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/
        return n_min * sqrt(Vd/300) * (mdot / 5e-6) * (1+ (n_max/n_min)*exp(-1.0*pow(((3.0*z/Lch) - 1.5),2)));
    };

    double Initial_Electron_Temperature(double z, double Te_Anode, double Te_Cathode, double Te_Max, double Lch, double z_max){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/
        return (1 - (z/z_max)) * Te_Anode + (z/z_max) * Te_Cathode + Te_Max * exp(-1.0*pow(((3.0*z/Lch) - 1.5),2));
    };

    double Initial_Ion_Bulk_Velocity(double Te_Anode, double Vd, double z, double Lch, double z_max){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/

        double u_Bohm = sqrt(1.6e-19 * Te_Anode / (1.6e-27 * 131.29));//e=1.6e-19 C, mi_Xe = 1.6e-27 (mp) * 131.29 (amu Xe)
        double u_max = sqrt(2 * 1.6e-19 * Vd / (1.6e-27 * 131.29));//e=1.6e-19 C, mi_Xe = 1.6e-27 (mp) * 131.29 (amu Xe)

        if (z<Lch){
            return u_Bohm + (2/3*(u_max - u_Bohm))*pow(z/Lch,2);
        } else{
            return (u_max + u_Bohm)*(1-(z-Lch)/(z_max-Lch)) + u_max * ((z-Lch)/(z_max-Lch));
        } 
    };


    HypiC::Particles_Object Initialize_Neutrals(HypiC::Options_Object Inputs)
    {
        int Particles_per_cell;
        double z;
        double dz;
        double ne;
        double Te;
        double v;
        double w;

        //calculate grid increment (assuming uniformly spaced)
        dz = Inputs.Domain_Length_m / Inputs.nCells;

        //create the particle class instance
        HypiC::Particles_Object Neutrals = HypiC::Particles_Object();
        Neutrals._IonizationDirection = -1.0;//neutrals are removed due to ionization

        //loop over grid cells
        for(size_t c=0; c<Inputs.nCells; ++c){
            //calculate number of particles per cell and the position
            Particles_per_cell = Inputs.N_Neutrals / Inputs.nCells;
            z = dz * c;
            //loop over the particles in the cell
            for(size_t p=0; p<Particles_per_cell; ++p){
                ne = HypiC::Initial_Electron_Density(z, Inputs.Initial_Min_Ion_Density,
                Inputs.Initial_Max_Ion_Density, Inputs.Discharge_Voltage_V, Inputs.Mass_Flow_Rate_kg_s, 
                Inputs.Channel_Length_m);
                Te = HypiC::Initial_Electron_Temperature(z, Inputs.Initial_Anode_Temperature_eV,
                Inputs.Initial_Cathode_Temperature_eV, Inputs.Initial_Max_Electron_Temperature_eV, 
                Inputs.Channel_Length_m, Inputs.Domain_Length_m);
                //calculate initial velocity
                //sample from maxwellian
                v=0;
                //calculate initial weight
                //see https://smileipic.github.io/Smilei/Understand/algorithms.html
                w = (ne / Particles_per_cell) * dz;

                //Add the particle
                Neutrals.Add_Particle(z, v, w, 0, ne, Te);
            }
        }

        //return
        return Neutrals;
    };

    HypiC::Particles_Object Initialize_Ions(HypiC::Options_Object Inputs)
    {
        int Particles_per_cell;
        double z;
        double dz;
        double ne;
        double Te;
        double v;
        double w;

        //calculate grid increment (assuming uniformly spaced)
        dz = Inputs.Domain_Length_m / Inputs.nCells;

        //create the particle class instance
        HypiC::Particles_Object Ions = HypiC::Particles_Object();
        Ions._IonizationDirection = 1.0;//ions are added due to ionization

        //loop over grid cells
        for(size_t c=0; c<Inputs.nCells; ++c){
            //calculate number of particles per cell and the position
            Particles_per_cell = Inputs.N_Ions / Inputs.nCells;
            z = dz * c;
            //loop over the particles in the cell
            for(size_t p=0; p<Particles_per_cell; ++p){
                ne = HypiC::Initial_Electron_Density(z, Inputs.Initial_Min_Ion_Density,
                Inputs.Initial_Max_Ion_Density, Inputs.Discharge_Voltage_V, Inputs.Mass_Flow_Rate_kg_s, 
                Inputs.Channel_Length_m);
                Te = HypiC::Initial_Electron_Temperature(z, Inputs.Initial_Anode_Temperature_eV,
                Inputs.Initial_Cathode_Temperature_eV, Inputs.Initial_Max_Electron_Temperature_eV, 
                Inputs.Channel_Length_m, Inputs.Domain_Length_m);
                //calculate initial velocity
                //sample from a maxwellian
                v=0;
                //calculate initial weight
                //see https://smileipic.github.io/Smilei/Understand/algorithms.html
                w = (ne / Particles_per_cell) * dz;

                //Add the particle
                Ions.Add_Particle(z, v, w, 0, ne, Te);
            }
        }

        //return
        return Ions;
    };

}