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
        int cell2;


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
            //actual pushing 
            Neutrals.Push_Particle(i, Simulation_Parameters.dt, Grid);

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

        //#pragma omp parallel for reduction(-:Neutrals._nParticles)
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

        //#pragma omp parallel for private(v,z,w)
        for (size_t i=0; i<n_remove; ++i){
            //sample from the maxwellian (force positive so resample)
            v = -1;
            while (v < 0){
                v = HypiC::Maxwellian_Sampler(un, u_Th);
            }
            //position is then sampled for distance from the anode 
            //see Birdsall pg 406
            z = ((double) rand()/RAND_MAX) * v * Simulation_Parameters.dt;

            //weights use the number denisty
            //might need to check the current cell n particles
            w = nn * (Simulation_Parameters.Domain_Length_m / Simulation_Parameters.nCells) / (Simulation_Parameters.N_Neutrals / Simulation_Parameters.nCells);

            //might need the grid methods for adding the electron quantities, however, this may get updated before the next
            //time step so not matter too much
            Neutrals.Add_Particle(z, v, w, 0, 1, 0);
        }

        n_remove = 0;//reset for ions

        //#pragma omp parallel for private(z) reduction(+:n_reflect)
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

        //#pragma omp parallel for private(current_index) reduction(+:Neutrals._nParticles)
        for (size_t i=n_reflect; i>0; --i){
            //add to the neutrals with position = 0 and reflected velocity
            //resample from Maxwellian?
            current_index = Reflect_These[i-1];
            Neutrals.Add_Particle(0, -1 * Ions.get_Position(current_index), Ions.get_Weight(current_index), 0, 1,
            Ions.get_ElectricField(current_index));
            Neutrals._nParticles += 1; //add to neutrals

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
       // #pragma omp parallel for private(current_index)
        for (size_t i=n_reflect; i>0; --i){
            //add to the neutrals with position = 0 and reflected velocity
            //resample from Maxwellian?
            current_index = Reflect_These[i-1];
            Neutrals.Add_Particle(0, -1 * Ions.get_Position(current_index), Ions.get_Weight(current_index), 0, 1, Ions.get_ElectricField(current_index));

            //remove the particle from the ions
            Ions.Remove_Particle(current_index);
            Remove_These.pop_back();//remove from the reflect list
        }

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