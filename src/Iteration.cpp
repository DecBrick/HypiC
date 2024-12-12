#include "HypiCpp.hpp"
#include <iostream>
namespace HypiC{
    HypiC::Particles_Object Update_Heavy_Species_Neutrals(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_rates, HypiC::Options_Object Simulation_Parameters){
        std::vector<size_t> Remove_These;
        std::vector<size_t> Reflect_These;
        size_t n_remove = 0;
        size_t n_reflect = 0;
        size_t current_index;
        double z;
        double v;
        double mass;
        double kb;
        double un;
        double u_Th;
        double nn;
        double w;


        //push the neutrals
        for (size_t i=0; i<Neutrals._nParticles; ++i){
            //push the neutrals
            Neutrals.Update_Particle(i, Simulation_Parameters.dt, Ionization_rates);

            //check outflow condition
            z = Neutrals.get_Position(i);
            if (z > Simulation_Parameters.Domain_Length_m){
                Remove_These.push_back(i);
            }
            //check and enforce reflection 
            if (z < 0){
                //set the position to 0 and switch sign on the velocity
                //resample from Maxwellian? 
                Neutrals.set_Position(i, 0.0);
                Neutrals.set_Velocity(i, -1.0 * Neutrals.get_Velocity(i)); 
            }
        }

        //enforce neutral boundary conditions
        for (size_t i=n_remove; i>0; --i){
            Neutrals.Remove_Particle(Remove_These[i]);//remove the particle
            Neutrals._nParticles -= 1;//update the number of particles
            Remove_These.pop_back();//remove from the remove list
        }
        
        //injection
        //sample from Maxwellian and place in the region using the netural velocity. 
        //First pass is to inject everytime a particle is removed
        mass = 131.29 * 1.66053907e-27;//for Xe
        kb = 1.380649e-23;
        un = sqrt(2 * kb * Simulation_Parameters.Initial_Neutral_Temperature_K / (M_PI * mass));
        u_Th = sqrt(kb * Simulation_Parameters.Initial_Neutral_Temperature_K / mass);
        nn = Simulation_Parameters.Mass_Flow_Rate_kg_s / (mass * un);
        srand(time(NULL));

        for (size_t i=0; i<n_remove; ++i){
            //sample from the maxwellian (force positive so resample)
            v = -1;
            while (v < 0){
                v = HypiC::Maxwellian_Sampler(un, u_Th);
            }
            //position is then sampled for distance from the anode 
            //see Birdsall pg 406
            z = (rand()/RAND_MAX) * v * Simulation_Parameters.dt;

            //weights use the number denisty
            //might need to check the current cell n particles
            w = nn * (Simulation_Parameters.Domain_Length_m / Simulation_Parameters.nCells) / (Simulation_Parameters.N_Neutrals / Simulation_Parameters.nCells);

            //might need the grid methods for adding the electron quantities, however, this may get updated before the next
            //time step so not matter too much
            Neutrals.Add_Particle(z, v, w, 0, 0, 0);
        }

        n_remove = 0;//reset for ions

        //push the ions
        for (size_t i=0; i<Ions._nParticles; ++i){


            //check outflow condition
            z = Ions.get_Position(i);

            if (z < 0){
                Reflect_These.push_back(i);
                n_reflect += 1;
            }
        }
        
        //enforce ion boundary conditions

        for (size_t i=n_reflect; i>0; --i){
            //add to the neutrals with position = 0 and reflected velocity
            //resample from Maxwellian?
            current_index = Reflect_These[i-1];
            Neutrals.Add_Particle(0, -1 * Ions.get_Position(current_index), Ions.get_Weight(current_index),
            Ions.get_ElectricField(current_index), Ions.get_ElectronDensity(current_index), Ions.get_ElectronTemperature(current_index));
            Neutrals._nParticles += 1; //add to neutrals

        }
        return Neutrals;
    };
    
    HypiC::Particles_Object Update_Heavy_Species_Ions(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_rates, HypiC::Options_Object Simulation_Parameters){
        std::vector<size_t> Remove_These;
        std::vector<size_t> Reflect_These;
        size_t n_remove = 0;
        size_t n_reflect = 0;
        size_t current_index;
        double z;
        double v;
        double mass;
        double kb;
        double un;
        double u_Th;
        double nn;
        double w;

        //push the ions
        for (size_t i=0; i<Ions._nParticles; ++i){
            Ions.Update_Particle(i, Simulation_Parameters.dt, Ionization_rates);

            //check outflow condition
            z = Ions.get_Position(i);
            if (z > Simulation_Parameters.Domain_Length_m){
                Remove_These.push_back(i);
                n_remove += 1;
            }
            if (z < 0){
                Reflect_These.push_back(i);
                n_reflect += 1;
            }
        }
        
        //enforce ion boundary conditions
        for (size_t i=n_remove; i>0; --i){
            Ions.Remove_Particle(Remove_These[i-1]);//remove the particle
            Ions._nParticles -= 1;//update the number of particles
            Remove_These.pop_back();//remove from the remove list
        }
        for (size_t i=n_reflect; i>0; --i){
            //add to the neutrals with position = 0 and reflected velocity
            //resample from Maxwellian?
            current_index = Reflect_These[i-1];
            //Neutrals.Add_Particle(0, -1 * Ions.get_Position(current_index), Ions.get_Weight(current_index),
            //Ions.get_ElectricField(current_index), Ions.get_ElectronDensity(current_index), Ions.get_ElectronTemperature(current_index));
            Neutrals._nParticles += 1; //add to neutrals

            //remove the particle from the ions
            Ions.Remove_Particle(current_index);
            Ions._nParticles -= 1;//update the number of particles
            Remove_These.pop_back();//remove from the reflect list
        }
        return Ions;
    };

    HypiC::Electrons_Object Update_Electrons(HypiC::Electrons_Object Electrons, HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_Rates, HypiC::Rate_Table_Object Collision_Loss_Rates, HypiC::Options_Object Simulation_Parameters){ //,double t){
        
        //first update the electron mobility and Temperature
        Electrons.Update_Mobility(Simulation_Parameters, Ionization_Rates);

        //then compute pressure gradient
        Electrons.Update_Pressure_Gradient(Simulation_Parameters);

        //Update discharge current
        Electrons.Compute_Discharge_Current(Simulation_Parameters);

        //Electron Velocity 
        Electrons.Update_Velocity(Simulation_Parameters);

        //Compute E Field, Potential
        Electrons.Compute_Electric_Field(Simulation_Parameters);
        Electrons.Solve_Potential(Simulation_Parameters);
        
        //update thermal conductivity and solve for energy
        Electrons.Update_Thermal_Conductivity(Simulation_Parameters);
        Electrons.Update_Electron_Energy(Simulation_Parameters, Collision_Loss_Rates);


        return Electrons;
    }

    

    
    //double Freq_Electron_Ion(double EnergyDensity, double ElectronTemp, double IonZ){
    //    return 2.9e-12 * pow(IonZ,2) * EnergyDensity * Coulomb_Logarithm(EnergyDensity, ElectronTemp, IonZ) / sqrt(pow(ElectronTemp,3));
    //}

    //double Coulomb_Logarithm(double EnergyDensity, double ElectronTemp, double IonZ){
    //    if (ElectronTemp < 10 * pow(IonZ,2)) {
    //        return 23 - 0.5 * log(1e-6 * EnergyDensity * pow(IonZ,2)/ pow(ElectronTemp,3));
    //    } else {
    //        return 24 - 0.5 * log(1e-6 * EnergyDensity / pow(ElectronTemp,2));
    //    }
    //}
}