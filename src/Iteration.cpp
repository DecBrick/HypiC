#include "HypiCpp.hpp"
#include <iostream>

namespace HypiC{
    HypiC::Particles_Object Update_Heavy_Species_Neutrals(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Grid, HypiC::Options_Object Simulation_Parameters){
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
        double t;
        double T;
        int cell2;
        double w_reinject = 0.0;
        double E_reinject = 0.0;
        double W_gen;
        double w_0_sum = 0.0;
        size_t N_0 = 0;
        size_t N_mass;
        size_t N_diffuse;
        size_t N_recombination;


        mass = 131.29 * 1.66053907e-27;//for Xe
        kb = 1.380649e-23;

        //#pragma omp parallel for private(z)
        //push the neutrals
        for (size_t i=0; i<Neutrals._nParticles; ++i){
            //first handle ionization for both cells the neutral is attached to 
            if (Ions._Ionization_Flag[Neutrals.get_Cell(i)]){
                //calculate particle contribution to the cell 
                nn = (1 - fabs((Neutrals.get_Position(i) - Grid.Get_CellCenter(Neutrals.get_Cell(i)))/Grid.Grid_Step)) 
                * Neutrals.get_Weight(i) / Grid.Grid_Step;
                /*if (Neutrals.get_Cell(i) == 125){
                    std::cout << "-----------\n";
                    std::cout << Neutrals.get_Position(i) << "\n";
                    std::cout << Grid.Get_CellCenter(Neutrals.get_Cell(i)) << "\n";
                    std::cout << Grid.Grid_Step << "\n";
                    std::cout << fabs((Neutrals.get_Position(i) - Grid.Get_CellCenter(Neutrals.get_Cell(i)))/Grid.Grid_Step) << "\n";
                    std::cout << nn << "\n";
                    std::cout << Grid.Get_NeutralDensity(Neutrals.get_Cell(i)) << "\n";
                    std::cout << Grid.Delta_ni[Neutrals.get_Cell(i)] << "\n";
                    std::cout <<  Neutrals.get_Weight(i) << "\n";
                }*/
                //weight is reduced by a proportional factor 
                w = Neutrals.get_Weight(i) - (nn/Grid.Get_NeutralDensity(Neutrals.get_Cell(i)))*Grid.Delta_ni[Neutrals.get_Cell(i)] / Grid.Grid_Step;
                Neutrals.set_Weight(i, w);
            }
            cell2 = Neutrals.get_Cell(i) + Neutrals.get_Cell2(i);
            if ((cell2 > 0) && (cell2 < Grid._nElectrons)){
                if (Ions._Ionization_Flag[cell2]){
                    //calculate particle contribution to the cell 
                    nn = (1 - fabs((Neutrals.get_Position(i) - Grid.Get_CellCenter(cell2)/Grid.Grid_Step))) 
                    * Neutrals.get_Weight(i) / Grid.Grid_Step;
                    //weight is reduced by a proportional factor 
                    w = Neutrals.get_Weight(i) - (nn/Grid.Get_NeutralDensity(cell2))*Grid.Delta_ni[cell2] / Grid.Grid_Step;
                    Neutrals.set_Weight(i, w);
                }
            }
            //total the neutrals in the first cell and time average
            //see Eq. 2.33 of Dominguez Vazquez thesis (2019)
            if ((Neutrals.get_Cell(i) == 0) || (cell2 == 0)){
                N_0 +=1;
                w_0_sum += Neutrals.get_Weight(i);
            }
            Neutrals.set_WBar((49 * Neutrals._WBar + w_0_sum) / 50);


            //actual pushing 
            Neutrals.Push_Particle(i, Simulation_Parameters.dt, Grid);

            //check outflow condition
            z = Neutrals.get_Position(i);
            if (z > Simulation_Parameters.Domain_Length_m){
                Remove_These.push_back(i);
            }
            //check and enforce reflection 
            if (z < 0){
                //assume completly diffuse reflection. Add weights and energy for this population
                w_reinject += Neutrals.get_Weight(i);
                E_reinject += Neutrals.get_Weight(i) * 0.5 * mass * pow(Neutrals.get_Velocity(i),2.0);

                //mark the current particle for removal 
                Remove_These.push_back(i);
            }
        }

        //#pragma omp parallel for reduction(-:Neutrals._nParticles)
        //enforce neutral boundary conditions
        for (size_t i=n_remove; i>0; --i){
            Neutrals.Remove_Particle(Remove_These[i]);//remove the particle
            Neutrals._nParticles -= 1;//update the number of particles
            Remove_These.pop_back();//remove from the remove list
        }
        
        //injection
        //have three populations: inflow, diffuse neutrals, and ion recombination
        //for the inflow, the area terms cancel so flux = mass flow/m 
        
        //first calculate generation weight
        W_gen = Neutrals._WBar * (N_0 / 100); //target 500 particles per cell for neutrals, 2x for 1D accounting (particle touches 2 cells)
        //now calculate number of particles generated for each population
        N_mass = ((size_t) (Simulation_Parameters.Mass_Flow_Rate_kg_s * Simulation_Parameters.dt / mass + Neutrals._wInject)/ (W_gen));
        N_diffuse = ((size_t) (w_reinject+Neutrals._wDiffuse) / W_gen);
        N_recombination = ((size_t) (Ions._ReinjectWeight+Neutrals._wRecombination) / W_gen);
        
        //tally remaining weights for next timestep 
        Neutrals.set_WInject((Simulation_Parameters.Mass_Flow_Rate_kg_s * Simulation_Parameters.dt / mass + Neutrals._wInject) - W_gen * N_mass);
        Neutrals.set_WDiffuse((w_reinject+Neutrals._wDiffuse) - W_gen * N_diffuse);
        Neutrals.set_WRecombination((Ions._ReinjectWeight+Neutrals._wRecombination) - W_gen * N_recombination);

        //sample for injection
        srand(time(NULL));
        T = Simulation_Parameters.Initial_Neutral_Temperature_K;//purly  wall temperature 
        for (size_t i=0; i<N_mass; ++i){
            //sample time fraction
            t = ((double)  rand()/RAND_MAX) * Simulation_Parameters.dt;
            //sample velocity, inject at thermal speed
            v = HypiC::Injection_Sampler(sqrt(2*kb * T / mass), T);
            //position follows
            z = v * t;
            //add the particle, with a switch for the side of the cell 
            if (z > 0.5 * Grid.Grid_Step){
                Neutrals.Add_Particle(z, v, W_gen, 0, 1, 0);
            }else{
                Neutrals.Add_Particle(z, v, W_gen, 0, -1, 0);
            }
        }

        //sample for diffusion
        T = E_reinject / 2.0;
        for (size_t i=0; i<N_diffuse; ++i){
            //sample time fraction
            t = ((double)  rand()/RAND_MAX) * Simulation_Parameters.dt;
            //sample velocity
            v = HypiC::Injection_Sampler(0, T);
            //position follows
            z = v * t;
            //add the particle, with a switch for the side of the cell 
            if (z > 0.5 * Grid.Grid_Step){
                Neutrals.Add_Particle(z, v, W_gen, 0, 1, 0);
            }else{
                Neutrals.Add_Particle(z, v, W_gen, 0, -1, 0);
            }
        }
        //sample for recombination
        T = Simulation_Parameters.Initial_Neutral_Temperature_K;//purly  wall temperature 
        for (size_t i=0; i<N_recombination; ++i){
            //sample time fraction
            t = ((double)  rand()/RAND_MAX) * Simulation_Parameters.dt;
            //sample velocity
            v = HypiC::Injection_Sampler(0, T);
            //position follows
            z = v * t;
            //add the particle, with a switch for the side of the cell 
            if (z > 0.5 * Grid.Grid_Step){
                Neutrals.Add_Particle(z, v, W_gen, 0, 1, 0);
            }else{
                Neutrals.Add_Particle(z, v, W_gen, 0, -1, 0);
            }
        }
        
        return Neutrals;
    };
    
    HypiC::Particles_Object Update_Heavy_Species_Ions(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Grid, HypiC::Options_Object Simulation_Parameters){
        std::vector<size_t> Remove_These;
        std::vector<size_t> Reflect_These;
        size_t n_remove = 0;
        size_t n_reflect = 0;
        size_t current_index;
        int c2;
        double z;
        double v;
        double mass;
        double kb;
        double un;
        double u_Th;
        double nn;
        double w;
        double z_rel;
        double s;
        double E_fld;
        double w_reinject;
        double E_reinject;


        mass = 131.29 * 1.66053907e-27;//for Xe

        srand(time(NULL));
        //ionization particle generation
        for (size_t i=0; i<Grid._nElectrons; ++i){
            //random chance for particle generation
            if ((Grid.Delta_ni[i]/Grid.Delta_ni_sum) > ((double) rand()/RAND_MAX)){
                //determine new particle properties
                z = Grid.Grid_Step * (i + ((double) rand()/RAND_MAX));
                v = HypiC::Maxwellian_Sampler(Grid.Neutral_Velocity_m_s[i], sqrt(kb * Grid.Neutral_Temperature_K[i] / mass));
                w = Grid.Delta_ni[i] / Grid.Grid_Step; 

                //check for cell2 index 
                if (z < Grid.Grid_Step*i + (Grid.Grid_Step / 2.0)){
                    c2 = -1;
                }else{
                    c2 = 1;
                }

                
                // shape factor (s for next cell = z_rel)
                z_rel = fabs((Grid.Cell_Center[i] - z)/Grid.Grid_Step);
                s = 1 - z_rel;
                // calculate partial electric field from this cell and next
                E_fld = s*Grid.Get_ElectricField(i);
                // add next cell contribution for all but last cell
                if((c2 > 0) && (c2 < Grid._nElectrons)){
                    E_fld += z_rel*Grid.Get_ElectricField(i+1);
                }

                //check for cell2 index 
                if (z < Grid.Grid_Step*i + (Grid.Grid_Step / 2.0)){
                    c2 = -1;
                }else{
                    c2 = 1;
                }

                //add new particle
                Ions.Add_Particle(z,v,w,i,c2, E_fld);

                //set flag for neutrals to be reduced
                Ions._Ionization_Flag[i] = true;

            }
            else{
                //set flag for neutrals to be reduced
                Ions._Ionization_Flag[i] = false;
            }
        }

        //#pragma omp parallel for private(z) reduction(+:n_remove,n_reflect)
        //push the ions
        for (size_t i=0; i<Ions._nParticles; ++i){
            Ions.Push_Particle(i, Simulation_Parameters.dt, Grid);
            if (Ions.get_Weight(i) <=0){
                std::cout << "Empty Particle\n";
            }
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
        //#pragma omp parallel for
        //count backwards on the remove so we don't mess up the indices 
        for (size_t i=n_remove; i>0; --i){
            Ions.Remove_Particle(Remove_These[i-1]);//remove the particle
            Remove_These.pop_back();//remove from the remove list
        }
        w_reinject = 0.0;
        E_reinject = 0.0;
       // #pragma omp parallel for private(current_index)
        for (size_t i=n_reflect; i>0; --i){
            //the particle to be reflected
            current_index = Reflect_These[i-1];
            //account for reinject weight and energy
            w_reinject += Ions.get_Weight(current_index);
            E_reinject += Ions.get_Weight(current_index) * (0.5*mass*pow(Ions.get_Velocity(current_index), 2.0) + 1.6e-19 * (Grid.Get_Potential(0) - Simulation_Parameters.Discharge_Voltage_V));


            //remove the particle from the ions
            Ions.Remove_Particle(current_index);
            Remove_These.pop_back();//remove from the reflect list
        }

        //assign the reinject quantities
        Ions.set_Reinjection(w_reinject, E_reinject / w_reinject);

        return Ions;
    };

    HypiC::Electrons_Object Update_Electrons(HypiC::Electrons_Object Electrons, HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Rate_Table_Object Ionization_Rates, HypiC::Rate_Table_Object Collision_Loss_Rates, HypiC::Options_Object Simulation_Parameters){ //,double t){
        
        //reset the ionization counter
        for (size_t i=0; i<Electrons._nElectrons; ++i){
            if (Ions._Ionization_Flag[i]){
                Electrons.Delta_ni[i] = 0.0;
            }
        }

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

}