#include <cmath>
#include <iostream>
#include "HypiCpp.hpp"

namespace HypiC{
    HypiC::Electrons_Object Particles_to_Grid(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons){
        // remove previous particle data to start summations at 0.0
        Electrons.Clear_Out_Particles();
        double dz = Electrons.Grid_Step;
        double z_p;
        double z_cell;
        double z_rel;
        double s; 
        double w; 
        double N_den;
        double N_flux;
        double N_energy;
        int cell1;
        int cell2;

        //#pragma omp parallel for //collapse(2)
        // loop over neutrals
        for(size_t i=0; i<Neutrals._nParticles; ++i){
            //pull the cells
            cell1 = Neutrals.get_Cell(i);
            cell2 = cell1 + Neutrals.get_Cell2(i);

            //shape factor using linear interpolation 
            z_p = Neutrals.get_Position(i);
            z_cell = Electrons.Get_CellCenter(cell1);
            z_rel = fabs( (z_cell - z_p)/dz);
            s = 1 - z_rel;
           
            // weight
            w = Neutrals.get_Weight(i);
            // calculate partial neutral density
            N_den = s*w;
            // calculate partial neutral flux for this cell and next
            N_flux = s*w*Neutrals.get_Velocity(i);
            // calculate energy term
            N_energy = s*w*pow(Neutrals.get_Velocity(i),2);
            // add values to cells' summations
            Electrons.Update_From_Neutrals(cell1,N_den,N_flux,N_energy);

            //update neighbor cell, avoid indexing errors   
            if((cell2 > 0) && (cell2 < Electrons._nElectrons)){
                // repeat for next cell
                // shape = z_rel
                N_den = z_rel*w;
                N_flux = z_rel*w*Neutrals.get_Velocity(i);
                N_energy = s*w*pow(Neutrals.get_Velocity(i),2);
                Electrons.Update_From_Neutrals(cell2,N_den,N_flux,N_energy);
            }
        }
        //#pragma omp parallel for //collapse(2)
        // loop over ions
        for(size_t i=0; i<Ions._nParticles; ++i){
            //pull the cells
            cell1 = Ions.get_Cell(i);
            cell2 = cell1 + Ions.get_Cell2(i);

            //shape factor using linear interpolation 
            z_p = Ions.get_Position(i);
            z_cell = Electrons.Get_CellCenter(cell1);
            z_rel = fabs( (z_cell - z_p)/dz);
            // shape factor (s for next cell = z_rel)
            s = 1 - z_rel;
            // weight
            w = Ions.get_Weight(i);
            // calculate partial plasma density for this cell and next
            N_den = s*w;
            // calculate partial current density for this cell and next
            N_flux = s*w*Ions.get_Velocity(i);
            // add values to cells' summations
            Electrons.Update_From_Ions(cell1,N_den,N_flux,N_flux);

            //update neighbor cell, avoid indexing errors   
            if((cell2 > 0) && (cell2 < Electrons._nElectrons)){
                // repeat for next cell
                // shape = z_rel
                N_den = z_rel*w;
                N_flux = z_rel*w*Ions.get_Velocity(i);
                Electrons.Update_From_Ions(cell2,N_den,N_flux,N_flux);
            }
        }
        //loop over cells to do the normalization using the total sums
        Electrons.Normalize_Interpolations();
        return Electrons;
    }
    /*//leaving this out for now, there is nothing the neutrals need from the electrons
    HypiC::Particles_Object Grid_to_Particles_Neutrals(HypiC::Particles_Object Neutrals, HypiC::Electrons_Object Electrons){
        #pragma omp parallel for //collapse(2)
        // loop over neutrals
        for(size_t i=0; i<Neutrals._nParticles; ++i){
            double e_den;
            double e_tmp;
            double dz = Electrons.Grid_Step;
            double z_p = Neutrals.get_Position(i);
            double z_cell;
            double z_rel; 
            double s;

            // loop over cells
            for(size_t c=0; c<Electrons._nElectrons; ++c){
                // determine if particle is within cell
                double z_cell = Electrons.Get_CellCenter(c);
                double z_rel = fabs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor (s for next cell = z_rel)
                    double s = 1 - z_rel;
                    // calculate partial electron density from this cell and next
                    e_den = s*Electrons.Get_PlasmaDensity(c);
                    // calculate partial electron temperature from this cell and next
                    e_tmp = s*Electrons.Get_ElectronTemperature(c);
                    // add next cell contribution for all but last cell
                    if(c!=Electrons._nElectrons){
                        e_den += z_rel*Electrons.Get_PlasmaDensity(c+1);
                        e_tmp += z_rel*Electrons.Get_ElectronTemperature(c+1);
                    }
                    break;
                }
            }

        }
        return Neutrals;
    }
    */

    HypiC::Particles_Object Grid_to_Particles_Ions(HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons){
        double E_fld;
        double dz = Electrons.Grid_Step;
        double z_p;
        double z_cell;
        double z_rel;
        double s;
        int cell;
        int cell2;
        
        #pragma omp parallel for //collapse(2)
        // loop over ions
        for(size_t i=0; i<Ions._nParticles; ++i){
            //grab cells
            cell = Ions.get_Cell(i);
            cell2 = cell + Ions.get_Cell2(i);
            //calculate shape factor
            z_p = Ions.get_Position(i);
            z_cell = Electrons.Get_CellCenter(cell);
            z_rel = fabs( (z_cell - z_p)/dz);
            s = 1 - z_rel;
            // calculate partial electric field from this cell and next
            E_fld = s*Electrons.Get_ElectricField(cell);
            // add next cell contribution for valid indices 
            //technically assumes the boundaries have no field 
            if((cell2 > 0) && (cell2 < Electrons._nElectrons)){
                E_fld += z_rel*Electrons.Get_ElectricField(cell + Ions.get_Cell2(i));
            }
            break;
            // set electric field
            Ions.set_ElectricField(i, E_fld);
        }
        return Ions;
    }
}
