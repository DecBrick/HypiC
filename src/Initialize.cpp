#include <cmath>
#include <fstream>
#include "HypiCpp.hpp"
#include <iostream>
#include <random>

namespace HypiC{
    double Maxwellian_Sampler(double mu, double sigma){
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution Dist{mu, sigma};
        //sample from a Maxwellian 
        return Dist(gen);
    };

    double Initial_Magnetic_Field(double B_max, double Lch, double z){
        double B;
        if (z < Lch){
            B = B_max * exp(-0.5 * pow(((z-Lch)/0.011),2.0));
        } else {
            B = B_max * exp(-0.5 * pow(((z-Lch)/0.018),2.0));
        }
        return B;
    }

    double Initial_Electron_Density(double z, double n_min, double n_max, double Vd, double mdot, double Lch){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/
        return  sqrt(Vd/300) * (mdot / 5e-6) * (n_min + (n_max - n_min)*exp(-1.0*pow(((3.0*z/Lch) - 1.5),2.0)));
    };

    double Initial_Neutral_Density(double z, double mass, double un, double mdot, double Lch, double Ach){
        double n_anode;
        double n_cathode;

        n_anode = mdot / (mass * un * Ach);
        n_cathode = 0.01 * n_anode;
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/
        return 0.5 * (n_anode + n_cathode + (n_anode - n_cathode) * tanh((z - 0.5 * Lch)/(Lch / 6.0)));
    };

    double Initial_Electron_Temperature(double z, double Te_Anode, double Te_Cathode, double Te_Max, double Lch, double z_max){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/explanation/initialization/
        double Te_min = std::min(Te_Anode, Te_Cathode);
        return (1 - (z/z_max)) * Te_Anode + (z/z_max) * Te_Cathode + (Te_Max - Te_min) * exp(-1.0*pow(((3.0*z/Lch) - 3.0),2.0));
    };

    double Initial_Ion_Bulk_Velocity(double Te_Anode, double Vd, double z, double Lch, double z_max){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/explanation/initialization/

        double u_Bohm = sqrt(1.6e-19 * Te_Anode / (1.6e-27 * 131.29));//e=1.6e-19 C, mi_Xe = 1.6e-27 (mp) * 131.29 (amu Xe)
        double u_max = sqrt(2 * 1.6e-19 * Vd / (1.6e-27 * 131.29));//e=1.6e-19 C, mi_Xe = 1.6e-27 (mp) * 131.29 (amu Xe)
        
        if (z<Lch){
            return u_Bohm + (2.0/3.0*(u_max - u_Bohm))*pow(z/Lch,2.0);
        } else{
            return (2* u_max / 3.0 + u_Bohm / 3.0)*(1-(z-Lch)/(z_max-Lch)) + u_max * ((z-Lch)/(z_max-Lch));
        } 
    };

    //for the initialization of particles, I think we can physically distribute the particles evenly
    //velocities come from maxwellian
    //weights come from density 
    //electron density and electron temperature come from the inital distributions
    //electric field can be set to 0 then updated (update electrons first?).

    HypiC::Particles_Object Initialize_Neutrals(HypiC::Options_Object Inputs)
    {
        int Particles_per_cell;
        double z;
        double z_particle;
        double dz;
        double nn;
        double v;
        double w;
        double mass;
        double kb;
        double un;
        int c2;

        //calculate grid increment (assuming uniformly spaced)
        dz = Inputs.Domain_Length_m / Inputs.nCells;

        //create the particle class instance
        HypiC::Particles_Object Neutrals = HypiC::Particles_Object();
        mass = 131.29 * 1.66053907e-27;//for Xe
        kb = 1.380649e-23;
        srand(time(NULL));

        #pragma omp parallel for //collapse(2)
        //loop over grid cells
        for(size_t c=0; c<Inputs.nCells; ++c){
            //calculate number of particles per cell and the position
            Particles_per_cell = Inputs.N_Neutrals / Inputs.nCells;
            z = dz * c;//left edge of cell 
            //loop over the particles in the cell
            for(size_t p=0; p<Particles_per_cell; ++p){
                //uniformly sample across cell. 
                z_particle = z + dz * (rand()/RAND_MAX);

                //calculate initial velocity
                //sample from maxwellian
                un = sqrt(2.0 * kb * Inputs.Initial_Neutral_Temperature_K / (M_PI * mass));
                v = HypiC::Maxwellian_Sampler(un,sqrt(kb * Inputs.Initial_Neutral_Temperature_K / mass));
                //calculate initial weight
                //see https://smileipic.github.io/Smilei/Understand/algorithms.html
                //use neutral density profile
                nn = Initial_Neutral_Density(z_particle, mass, un, Inputs.Mass_Flow_Rate_kg_s, Inputs.Channel_Length_m,
                Inputs.Channel_Area_m2);
                w = (nn / Particles_per_cell) * dz;
                
                //check for cell2 index 
                if (z_particle < z + (dz / 2.0)){
                    c2 = -1;
                }else{
                    c2 = 1;
                }
                //Add the particle
                Neutrals.Add_Particle(z_particle, v, w, c, c2, 0.0);
            }
        }

        //return
        return Neutrals;
    };

    HypiC::Particles_Object Initialize_Ions(HypiC::Options_Object Inputs)
    {
        int Particles_per_cell;
        double z;
        double z_particle;
        double dz;
        double ne;
        double v;
        double w;
        double mass;
        double kb;
        int c2;

        //calculate grid increment (assuming uniformly spaced)
        dz = Inputs.Domain_Length_m / Inputs.nCells;

        //create the particle class instance
        HypiC::Particles_Object Ions = HypiC::Particles_Object();
        mass = 131.29 * 1.66053907e-27;//for Xe
        kb = 1.380649e-23;
        Ions._ChargetoMassRatio = 1.602176634e-19 / mass; 
        srand(time(NULL));

        #pragma omp parallel for private(z,z_particle,ne,Te,v,w) //collapse(2)
        //loop over grid cells
        for(size_t c=0; c<Inputs.nCells; ++c){
            //calculate number of particles per cell and the position
            Particles_per_cell = Inputs.N_Ions / Inputs.nCells;
            z = dz * c;
            
            //loop over the particles in the cell
            for(size_t p=0; p<Particles_per_cell; ++p){
                //uniformly sample
                //later add to the position 
                z_particle = z + dz * (rand()/RAND_MAX);

                ne = HypiC::Initial_Electron_Density(z_particle, Inputs.Initial_Min_Ion_Density,
                Inputs.Initial_Max_Ion_Density, Inputs.Discharge_Voltage_V, Inputs.Mass_Flow_Rate_kg_s, 
                Inputs.Channel_Length_m);
                
                //calculate initial velocity
                //sample from a maxwellian
                v = HypiC::Maxwellian_Sampler(HypiC::Initial_Ion_Bulk_Velocity(Inputs.Initial_Anode_Temperature_eV, Inputs.Discharge_Voltage_V,
                z_particle, Inputs.Channel_Length_m, Inputs.Domain_Length_m), sqrt(kb * Inputs.Initial_Ion_Temperature_K / mass));
                if (v>3e8){
                    std::cout << z_particle << "\n";
                    std::cout << HypiC::Initial_Ion_Bulk_Velocity(Inputs.Initial_Anode_Temperature_eV, Inputs.Discharge_Voltage_V,
                    z_particle, Inputs.Channel_Length_m, Inputs.Domain_Length_m) << "\n";
                    
                }
                //calculate initial weight, for singly charged only can use the electron number density
                //see https://smileipic.github.io/Smilei/Understand/algorithms.html
                w = (ne / Particles_per_cell) * dz;

                //check for cell2 index 
                if (z_particle < z + (dz / 2.0)){
                    c2 = -1;
                }else{
                    c2 = 1;
                }

                //Add the particle
                Ions.Add_Particle(z, v, w, c, c2, 0);
            }

            //add an ionization flag
            Ions._Ionization_Flag.push_back(false);
        }

        //return
        return Ions;
    };

    HypiC::Electrons_Object Initialize_Electrons(HypiC::Options_Object Inputs)
    {
        int Particles_per_cell;
        double z;
        double z_particle;
        double dz;
        double ne;
        double Te;
        double v;
        double w;
        double mass;
        double kb;
        double EnergyDensity;
        double B;
        double f;
        double Efield;

        //calculate grid increment (assuming uniformly spaced)
        dz = Inputs.Domain_Length_m / Inputs.nCells;

        //create the particle class instance
        HypiC::Electrons_Object Electrons = HypiC::Electrons_Object();
        Electrons.Grid_Step = dz;
        mass = 131.29 * 1.66053907e-27;//for Xe
        kb = 1.380649e-23;
        srand(time(NULL));

        #pragma omp parallel for private(z,ne,Te,v,EnergyDensity,B,f,Efield)
        //loop over grid cells
        for(size_t c=0; c<Inputs.nCells; ++c){
            //calculate number of particles per cell and the position
            z = 0.5 * dz + dz * c;

            ne = HypiC::Initial_Electron_Density(z, Inputs.Initial_Min_Ion_Density,
            Inputs.Initial_Max_Ion_Density, Inputs.Discharge_Voltage_V, Inputs.Mass_Flow_Rate_kg_s, 
            Inputs.Channel_Length_m);
            Te = HypiC::Initial_Electron_Temperature(z, Inputs.Initial_Anode_Temperature_eV,
            Inputs.Initial_Cathode_Temperature_eV, Inputs.Initial_Max_Electron_Temperature_eV, 
            Inputs.Channel_Length_m, Inputs.Domain_Length_m);
            B = HypiC::Initial_Magnetic_Field(0.015, Inputs.Channel_Length_m, z);
            EnergyDensity = ne * Te * 3.0/2.0;

            v = HypiC::Initial_Ion_Bulk_Velocity(Inputs.Initial_Anode_Temperature_eV, Inputs.Discharge_Voltage_V,
            z, Inputs.Channel_Length_m, Inputs.Domain_Length_m);

            // Utilizing Case 2 from 
            // https://0534de96-08f4-4c2b-82e3-b2a3f3216551.filesusr.com/ugd/8243e7_030af6befce7412ca1232089b304903e.pdf
            if (z <= Inputs.Channel_Length_m) {
                f = 0.5e7;
            } else {
                f = 1e7;
            }

            Efield = 2;

            //Add the particle
            Electrons.Add_Electron(ne, Te, B, EnergyDensity, v, f, Efield, z);
        }

        //return
        return Electrons;
    };

}