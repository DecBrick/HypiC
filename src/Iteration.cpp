#include "HypiCpp.hpp"

namespace HypiC{
    void Update_Heavy_Species(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_rates, HypiC::Options_Object Simulation_Parameters){
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
            Ions.Remove_Particle(Remove_These[i]);//remove the particle
            Ions._nParticles -= 1;//update the number of particles
            Remove_These.pop_back();//remove from the remove list
        }
        for (size_t i=n_reflect; i>0; --i){
            //add to the neutrals with position = 0 and reflected velocity
            //resample from Maxwellian?
            current_index = Reflect_These[i];
            Neutrals.Add_Particle(0, -1 * Ions.get_Position(current_index), Ions.get_Weight(current_index),
            Ions.get_ElectricField(current_index), Ions.get_ElectronDensity(current_index), Ions.get_ElectronTemperature(current_index));
            Neutrals._nParticles += 1; //add to neutrals

            //remove the particle from the ions
            Ions.Remove_Particle(current_index);
            Ions._nParticles -= 1;//update the number of particles
            Remove_These.pop_back();//remove from the reflect list
        }
    };



    void Update_Electrons(HypiC::Electrons_Object Electrons, HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_Rates, HypiC::Options_Object Simulation_Parameters){ //,double t){
        for(size_t c=0; c<Simulation_Parameters.nCells; ++c){
            Electrons.Electron_Temperature_eV[c] = 2/3 * Electrons.EnergyDensity[c]/Electrons.Plasma_Density_m3[c];
            Electrons.Electron_Pressure[c] =  1.5 * Electrons.EnergyDensity[c] * Electrons.Electron_Temperature_eV[c];

            //Update Electron-Ion Collisions - NOT NECESSARY FOR LANDMARK AAAAA
            //Electrons.Freq_Elec_Ion[c] = Freq_Electron_Ion(Electrons.EnergyDensity[c],Electrons.Electron_Temperature_eV[c],Electrons.Ion_Z[c]);

            //Update OTHER
            //This seems to be super complicated?? But Landmark seems to be constant 2.5e-13 rate?
            Electrons.Freq_Elec_Neutral[c] = 2.5e-13 * Electrons.Neutral_Density_m3[c];
            
            //Electrons.Freq_Classical[c] = Electrons.Freq_Elec_Neutral[c] + Electrons.Freq_Elec_Ion[c];
            Electrons.Freq_Classical[c] = Electrons.Freq_Elec_Neutral[c];

            //Seems like Radial Loss Freq is a constant 1e7 for Landmark?
          
            Electrons.Freq_Electron_Wall_Collision[c] = 1e7 * Linear_Transition(Electrons.Cell_Center[c],Simulation_Parameters.Channel_Length_m, 0.2*Simulation_Parameters.Channel_Length_m,1,0);
        }
            //Anomalous transport?
        for(size_t c1=0; c1<Simulation_Parameters.nCells; ++c1){

            double Elec_Cycl_Freq = 1.602176634e-19 * Electrons.Magnetic_Field_G[c1] / 9.10938356e-31;

            //Might need to change this to 0.1, 1?
            double Beta = Linear_Transition(Electrons.Cell_Center[c1],Simulation_Parameters.Channel_Length_m, 0.2*Simulation_Parameters.Channel_Length_m, 1/160,1/16);

            Electrons.Freq_Anomalous_Collision[c1] = Elec_Cycl_Freq * Beta;

        }
            
        //Smooth? - Not necessary
        //std::vector<double> Freq_Anomalous_Copy = Electrons.Freq_Anomalous_Collision;
        //Smooth(Electrons.Freq_Anomalous_Collision, Freq_Anomalous_Copy, Iterations)
        for(size_t c1=0; c1<Simulation_Parameters.nCells; ++c1){
            Electrons.Freq_Total_Electron_Collision[c1] = Electrons.Freq_Electron_Wall_Collision[c1] + Electrons.Freq_Anomalous_Collision[c1] + Electrons.Freq_Classical[c1];

            double Omega =  1.602176634e-19 * Electrons.Magnetic_Field_G[c1] / (9.10938356e-31 *Electrons.Freq_Total_Electron_Collision[c1]);
            Electrons.Electron_Mobility[c1] = 1.602176634e-19 / (9.10938356e-31 * Electrons.Freq_Total_Electron_Collision[c1] * (1+pow(Omega,2)));
        
        }
        for (size_t c=0; c<Simulation_Parameters.nCells; ++c){
            Electrons.Ion_Current_Density[c] = Electrons.Ion_Z[c] * Electrons.Plasma_Density_m3[c] * Electrons.Ion_Velocity_m_s[c];
        }
        //Discharge Voltage
        double Discharge_Current = Integrate_Discharge_Current(Electrons, Simulation_Parameters);

        //Electron Velocity + Electron Kinetic
        for(size_t c1=0; c1<Simulation_Parameters.nCells; ++c1){
            Electrons.Electron_Velocity_m_s[c1] = (Electrons.Ion_Current_Density[c1] - Discharge_Current/Simulation_Parameters.Channel_Area_m2/1.602176634e-19 / Electrons.EnergyDensity[c1]);
            Electrons.Electron_Kinetic_Energy[c1] = 0.5 * 9.10938356e-31 * (1 + pow(1.602176634e-19*Electrons.Magnetic_Field_G[c1]/9.10938356e-31/Electrons.Freq_Total_Electron_Collision[c1],2)*pow(Electrons.Electron_Velocity_m_s[c1],2)/1.602176634e-19);
        }

        //Compute Pressure graidient, EFIeld, Potential, Thermal Conductivity, Energy
        Compute_Pressure_Gradient(Electrons,Simulation_Parameters);

        Compute_Electric_Field(Electrons,Simulation_Parameters, Discharge_Current);

        //Solve_Potential(Electrons,Simulation_Parameters);

        //Update_Thermal_Conductivity(Electrons,Simulation_Parameters);

        //Update_Electron_Energy(Electrons,Simulation_Parameters);
    }



    void Compute_Electric_Field(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters, double Discharge_Current){
        for(size_t i=0; i<Simulation_Parameters.nCells; ++i){
            double E = ((Discharge_Current / Simulation_Parameters.Channel_Area_m2 - Electrons.Ion_Current_Density[i])/1.602176634e-19/Electrons.Electron_Mobility[i] - Electrons.Electron_Pressure_Gradient[i]/Electrons.EnergyDensity[i]);
            Electrons.Electric_Field_V_m[i] = -E;
        }
    }

    void Compute_Pressure_Gradient(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters){
        Electrons.Electron_Pressure_Gradient[0] = Forward_Difference(Electrons.Electron_Pressure[0],Electrons.Electron_Pressure[1],Electrons.Electron_Pressure[2],Electrons.Cell_Center[0],Electrons.Cell_Center[1],Electrons.Cell_Center[2]);
        for(size_t i=1; i<Simulation_Parameters.nCells-1; ++i){
            Electrons.Electron_Pressure_Gradient[i] = Central_Difference(Electrons.Electron_Pressure[i-1],Electrons.Electron_Pressure[i],Electrons.Electron_Pressure[i+1],Electrons.Cell_Center[i-1],Electrons.Cell_Center[i],Electrons.Cell_Center[i+1]);
        }
        size_t end = Electrons.Electron_Pressure_Gradient.size() - 1;
        Electrons.Electron_Pressure_Gradient[end] = Backward_Difference(Electrons.Electron_Pressure[end-2],Electrons.Electron_Pressure[end-1],Electrons.Electron_Pressure[end],Electrons.Cell_Center[end-2],Electrons.Cell_Center[end-1],Electrons.Cell_Center[end]);
    }

    double Forward_Difference(double f0, double f1, double f2, double x0, double x1, double x2){
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = -(2*h1+h2)/h1/(h1+h2);
        double c1 = (h1+h2)/(h1*h2);
        double c2 = -h1/h2/(h1+h2);
        return c0 * f0 + c1 * f1 + c2 * f2;
    }

    double Central_Difference(double f0, double f1, double f2, double x0, double x1, double x2){
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = -h2/h1/(h1+h2);
        double c1 = -(h1+h2)/(h1*h2);
        double c2 = h1/h2/(h1+h2);
        return c0 * f0 + c1 * f1 + c2 * f2;
    }

    double Backward_Difference(double f0, double f1, double f2, double x0, double x1, double x2){
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = h2/h1/(h1+h2);
        double c1 = -(h1+h2)/(h1*h2);
        double c2 = (h1+2*h2)/h2/(h1+h2);

        return c0 * f0 + c1 * f1 + c2 * f2;
    }

    double Integrate_Discharge_Current(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters){
        double int1 = 0;
        double int2 = 0;
        double Dz = Simulation_Parameters.Domain_Length_m / Simulation_Parameters.nCells;

        for (size_t c=0; c<Simulation_Parameters.nCells-1; ++c){
            double dz = Dz * c;
            double int1_1 = (Electrons.Ion_Current_Density[c]/1.602176634e-19/Electrons.Electron_Mobility[c] + Electrons.Electron_Pressure_Gradient[c]/Electrons.EnergyDensity[c]);
            double int1_2 = (Electrons.Ion_Current_Density[c+1]/1.602176634e-19/Electrons.Electron_Mobility[c+1] + Electrons.Electron_Pressure_Gradient[c+1]/Electrons.EnergyDensity[c+1]);

            int1 += 0.5 * dz * (int1_1 + int1_2);

            double int2_1 = 1/(1.602176634e-19 * Electrons.EnergyDensity[c] * Electrons.Electron_Mobility[c]*Simulation_Parameters.Channel_Area_m2);
            double int2_2 = 1/(1.602176634e-19 * Electrons.EnergyDensity[c+1] * Electrons.Electron_Mobility[c+1]*Simulation_Parameters.Channel_Area_m2);

            int2 += 0.5 * dz * (int2_1 + int2_2);

        }
        double Discharge_Current = (Simulation_Parameters.Discharge_Voltage_V + int1)/int2;
        return Discharge_Current;
    }

    double Linear_Transition(double x, double cutoff, double L, double y1, double y2){
        double x1 = cutoff - L/2;
        double x2 = cutoff + L/2;
        if (x < x1){
            return y1;
        } else if (x > x2){
            return y2;
        } else {
            double t = (x-x1)/(x2-x1);
            return t * (y2-y1) + y1;
        }
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