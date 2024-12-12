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

        //Discharge Voltage
        double Discharge_Current = Integrate_Discharge_Current(Electrons, Simulation_Parameters);
        //std::cout << Discharge_Current << "\n";
        //Electron Velocity + Electron Kinetic
        for(size_t c1=0; c1<Simulation_Parameters.nCells; ++c1){
            Electrons.Electron_Velocity_m_s[c1] = (Electrons.Ion_Current_Density[c1] - Discharge_Current/Simulation_Parameters.Channel_Area_m2/1.602176634e-19 / Electrons.EnergyDensity[c1]);
            Electrons.Electron_Kinetic_Energy[c1] = 0.5 * 9.10938356e-31 * (1 + pow(1.602176634e-19*Electrons.Magnetic_Field_T[c1]/9.10938356e-31/Electrons.Freq_Total_Electron_Collision[c1],2)*pow(Electrons.Electron_Velocity_m_s[c1],2)/1.602176634e-19);
        }

        //Compute E Field, Potential, Thermal Conductivity, Energy
        Compute_Electric_Field(Electrons,Simulation_Parameters, Discharge_Current);

        Solve_Potential(Electrons,Simulation_Parameters);
        
        Update_Thermal_Conductivity(Electrons,Simulation_Parameters);

        Update_Electron_Energy(Electrons,Simulation_Parameters, Collision_Loss_Rates);
        return Electrons;
    }

    void Update_Electron_Energy(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters, HypiC::Rate_Table_Object Loss_Rates){
        std::vector<double> diag(Simulation_Parameters.nCells,1);
        std::vector<double> diag_low(Simulation_Parameters.nCells-1,1);
        std::vector<double> diag_up(Simulation_Parameters.nCells-1,1);
        std::vector<double> B(Simulation_Parameters.nCells,1);
        std::vector<double> Energy_new;
        double je_sheath;
        double T0;

        //setting boundary conditions
        diag[0] = 1;
        diag_up[0] = 0;
        diag[Simulation_Parameters.nCells-1] = 1;
        diag_low[Simulation_Parameters.nCells-2] = 0;

        //anode neuman for sheath boundary        
        B[0] = 0;
        diag[0] = 1/Electrons.EnergyDensity[0];
        diag_up[0] = -1/Electrons.EnergyDensity[1];

        //cathode boundary condition
        B[Simulation_Parameters.nCells-1] = 1.5 * Simulation_Parameters.Initial_Cathode_Temperature_eV * Electrons.EnergyDensity[Simulation_Parameters.nCells-1];

        //calculate source terms
        //there should be a timestep term in here too 
        std::vector<double> OhmicHeating(Simulation_Parameters.nCells,0);
        std::vector<double> Collisional_Loss(Simulation_Parameters.nCells,0);
        std::vector<double> WallPowerLoss(Simulation_Parameters.nCells,0);
        for(size_t i=0; i<Simulation_Parameters.nCells-1; ++i){ //-1 is b/c we want to fix the Cathode cell using a dirichlet condition 
            OhmicHeating[i] = -1.602176634e-19 * Electrons.Plasma_Density_m3[i] * Electrons.Electron_Velocity_m_s[i] * Electrons.Electric_Field_V_m[i];
            Collisional_Loss[i] = Electrons.Plasma_Density_m3[i] * Electrons.Neutral_Density_m3[i] * Loss_Rates.interpolate(Electrons.Electron_Temperature_eV[i]);
            
            if (Electrons.Cell_Center[i] <= Simulation_Parameters.Channel_Length_m){
                WallPowerLoss[i] = 7.5e6 * Electrons.Electron_Temperature_eV[i] * exp( - 40 / (3 *Electrons.Electron_Temperature_eV[i]));
            } else{
                WallPowerLoss[i] = 0;
            }
            B[i] = OhmicHeating[i] - Electrons.EnergyDensity[i] * WallPowerLoss[i] - Collisional_Loss[i] + Simulation_Parameters.dt * Electrons.EnergyDensity[i]; 
        }

        //calculate fluxes (setting up linear matrix)
        /*we have two terms we need to care about here are the convective flux (div 5/3 * flux*energy)
        and the heat flux (div kappa nabla energy). By applying the divergence theorem we can remove the div to 
        make the equation only based on the fluxes at the edges

        For the convective flux, use upwinding, which means that we use the value in the same direction as the velocity
        for 1D, if velocity is positive, use the cell to the left, if negative, use the cell to the right
        do this for both edges of each cell. Then just need a factor of 5/3*ue for the respective cell's energy density

        For the heat flux, use upwinding for which cell's kappa to apply and then use a central difference scheme to finite
        difference the gradient in the energy density which adds a geometric factor. 
        */

        //loop over main body
        for(size_t i=1; i<Simulation_Parameters.nCells - 1; ++i){
            if (Electrons.Electron_Velocity_m_s[i] > 0){ 
                //upwind from 0 to 1 and 1 to 2
                diag_low[i-1] = (5.0/3.0) * Electrons.Electron_Velocity_m_s[i-1] + Electrons.Electron_Thermal_Conductivity[i-1] / Electrons.Grid_Step;
                diag[i] = (5.0/3.0) * Electrons.Electron_Velocity_m_s[i] - (Electrons.Electron_Thermal_Conductivity[i-1] - Electrons.Electron_Thermal_Conductivity[i]) / Electrons.Grid_Step;
                diag_up[i] = -Electrons.Electron_Thermal_Conductivity[i] / Electrons.Grid_Step;
            } else{
                //upwind from 1 to 0 and 2 to 1
                diag_low[i-1] =  Electrons.Electron_Thermal_Conductivity[i] / Electrons.Grid_Step;
                diag[i] = (5.0/3.0) * Electrons.Electron_Velocity_m_s[i] - (Electrons.Electron_Thermal_Conductivity[i] - Electrons.Electron_Thermal_Conductivity[i+1]) / Electrons.Grid_Step;
                diag_up[i] = (5.0/3.0) * Electrons.Electron_Velocity_m_s[i+1] -Electrons.Electron_Thermal_Conductivity[i+1] / Electrons.Grid_Step;
            }

        }

        //apply sheath boundary condition 
        //neglect heat flux at left edge
        if (Electrons.Electron_Velocity_m_s[0] > 0){ 
            //upwind from 1 to 2
            diag[0] = (5.0/3.0) * Electrons.Electron_Velocity_m_s[0] + Electrons.Electron_Thermal_Conductivity[0] / Electrons.Grid_Step;
            diag_up[0] = -Electrons.Electron_Thermal_Conductivity[0] / Electrons.Grid_Step;
        } else{
            //upwind from 2 to 1
            diag[0] = Electrons.Electron_Thermal_Conductivity[1] / Electrons.Grid_Step;
            diag_up[0] = (5.0/3.0) * Electrons.Electron_Velocity_m_s[1] - Electrons.Electron_Thermal_Conductivity[1] / Electrons.Grid_Step;
        }

        //handle the left edge
        T0 = (2.0/3.0) * Electrons.EnergyDensity[0] / Electrons.Plasma_Density_m3[0];
        je_sheath = (Integrate_Discharge_Current(Electrons, Simulation_Parameters) / Simulation_Parameters.Channel_Area_m2) - Electrons.Ion_Current_Density[0];
        //the back term (1-log) of this is the sheath potential 
        diag[0] = (4.0/3.0) * je_sheath / (-1.602176634e-19 * Electrons.Plasma_Density_m3[0]) * (1.0 - log(std::min(1.0, je_sheath / (1.602176634e-19 * Electrons.Plasma_Density_m3[0] * sqrt(8 * 1.602176634e-19 * T0/ (M_PI * 9.10938356e-31)) / 4))));

        //call matrix solver, might want to use Thomas https://www.quantstart.com/articles/Tridiagonal-Matrix-Algorithm-Thomas-Algorithm-in-C/
        Energy_new = HypiC::Thomas_Algorithm(diag_low, diag, diag_up, B);

        //limit to a minimum temperature and assign
        for(size_t i=0; i<Simulation_Parameters.nCells; ++i){
            Electrons.EnergyDensity[i] = std::max(Energy_new[i], 1.5 * Electrons.Plasma_Density_m3[i] * Simulation_Parameters.Min_Electron_Temperature_eV);
        }
    }

    void Solve_Potential(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters){
        Electrons.Potential[0] = Simulation_Parameters.Discharge_Voltage_V;

        for(size_t i=1; i<Simulation_Parameters.nCells; ++i){
            double dx = Electrons.Cell_Center[i] - Electrons.Cell_Center[i-1];
            Electrons.Potential[i] = Electrons.Potential[i-1] + 0.5 * dx * (Electrons.Electric_Field_V_m[i] + Electrons.Electric_Field_V_m[i-1]);
        }
    }

    void Compute_Electric_Field(HypiC::Electrons_Object Electrons,HypiC::Options_Object Simulation_Parameters, double Discharge_Current){
        for(size_t i=0; i<Simulation_Parameters.nCells; ++i){
            double E = ((Discharge_Current / Simulation_Parameters.Channel_Area_m2 - Electrons.Ion_Current_Density[i])/1.602176634e-19/Electrons.Electron_Mobility[i] - Electrons.Electron_Pressure_Gradient[i]/Electrons.EnergyDensity[i]);
            Electrons.Electric_Field_V_m[i] = -E;
        }
    }

    void Update_Thermal_Conductivity(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters){
        for(size_t i=0; i<Simulation_Parameters.nCells; ++i){
            Electrons.Electron_Thermal_Conductivity[i] = (10.0 / (9.0 * 1.602176634e-19)) * Electrons.Electron_Mobility[i] * Electrons.Plasma_Density_m3[i] * Electrons.EnergyDensity[i];
        }
    }


    double Integrate_Discharge_Current(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters){
        double int1 = 0;
        double int2 = 0;
        double Dz = Simulation_Parameters.Domain_Length_m / Simulation_Parameters.nCells;
        double dz;
        double int1_1;
        double int1_2;
        double int2_1;
        double int2_2;
        double enemu1;
        double enemu2;

        for (size_t c=0; c<Simulation_Parameters.nCells-1; ++c){
            //dz = Dz * c;
            enemu1 = 1.602176634e-19 * Electrons.Plasma_Density_m3[c] * Electrons.Electron_Mobility[c];
            enemu2 = 1.602176634e-19 * Electrons.Plasma_Density_m3[c+1] * Electrons.Electron_Mobility[c+1];
            int1_1 = (Electrons.Ion_Current_Density[c]/enemu1 + Electrons.Electron_Pressure_Gradient[c]/Electrons.Plasma_Density_m3[c]/1.602176634e-19 );
            int1_2 = (Electrons.Ion_Current_Density[c+1]/enemu2 + Electrons.Electron_Pressure_Gradient[c+1]/Electrons.Plasma_Density_m3[c+1]/1.602176634e-19);

            int1 += 0.5 * Dz * (int1_1 + int1_2);

            //the area is only correct within the channel, need to add plume area model
            int2_1 = 1/(enemu1*Simulation_Parameters.Channel_Area_m2);
            int2_2 = 1/(enemu2*Simulation_Parameters.Channel_Area_m2);

            int2 += 0.5 * Dz * (int2_1 + int2_2);

        }

        double Discharge_Current = (Simulation_Parameters.Discharge_Voltage_V + int1)/int2;
        return Discharge_Current;
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