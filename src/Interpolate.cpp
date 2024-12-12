#include <cmath>
#include <iostream>
#include "HypiCpp.hpp"

namespace HypiC{
    HypiC::Electrons_Object Particles_to_Grid(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons){
        // remove previous particle data to start summations at 0.0
        Electrons.Clear_Out_Particles(Neutrals._nParticles,Ions._nParticles);
        
        double dz = Electrons.Grid_Step;
        // loop over neutrals
        for(size_t i=0; i<Neutrals._nParticles; ++i){
            double z_p = Neutrals.get_Position(i);

            // loop over cells
            for(size_t c=0; c<Electrons._nElectrons; ++c){
                // determine if particle is within cell
                double z_cell = Electrons.Get_CellCenter(c);
                double z_rel = fabs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor (s for next cell = z_rel)
                    double s = 1 - z_rel;
                    // weight
                    double w = Neutrals.get_Weight(i);
                    // calculate partial neutral density for this cell and next
                    double N_den = s*w /dz;
                    double N_den_next = z_rel*w/dz;
                    // calculate partial neutral flux for this cell and next
                    double N_flux = s*(w/dz)*Neutrals.get_Velocity(i);
                    double N_flux_next = z_rel*(w/dz)*Neutrals.get_Velocity(i);
                    // add values to cells' summations
                    Electrons.Update_From_Neutrals(c,N_den,N_den_next,N_flux/N_den,N_flux_next/N_den_next);
                    break;
                }
            }
        }

        // loop over ions
        for(size_t i=0; i<Ions._nParticles; ++i){
            double z_p = Ions.get_Position(i);

            // loop over cells
            for(size_t c=0; c<Electrons._nElectrons; ++c){
                // determine if particle is within cell
                double z_cell = Electrons.Get_CellCenter(c);
                double z_rel = fabs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor (s for next cell = z_rel)
                    double s = 1 - z_rel;
                    // weight
                    double w = Ions.get_Weight(i);
                    // calculate partial plasma density for this cell and next
                    double P_den = s*w /dz;
                    double P_den_next = z_rel*w/dz;
                    // calculate partial current density for this cell and next
                    double J_den = s*(w/dz)*Ions.get_Velocity(i);
                    double J_den_next = z_rel*(w/dz)*Ions.get_Velocity(i);
                    // add values to cells' summations
                    Electrons.Update_From_Ions(c,P_den,P_den_next,1.602176634e-19*J_den,1.602176634e-19*J_den_next,J_den/P_den,J_den_next/P_den_next);
                    break;
                }
            }
        }
        return Electrons;
    }

    HypiC::Particles_Object Grid_to_Particles_Neutrals(HypiC::Particles_Object Neutrals, HypiC::Electrons_Object Electrons){
        // loop over neutrals
        for(size_t i=0; i<Neutrals._nParticles; ++i){
            double e_den;
            double e_tmp;
            double dz = Electrons.Grid_Step;
            double z_p = Neutrals.get_Position(i);

            // loop over cells
            for(size_t c=0; c<Electrons._nElectrons; ++c){
                // determine if particle is within cell
                double z_cell = Electrons.Get_CellCenter(c);
                double z_rel = abs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor (s for next cell = z_rel)
                    double s = 1 - z_rel;
                    // calculate partial electron density from this cell and next
                    e_den = s*Electrons.Get_PlasmaDensity(c) + z_rel*Electrons.Get_PlasmaDensity(c+1);
                    // calculate partial electron temperature from this cell and next
                    e_tmp = s*Electrons.Get_ElectronTemperature(c) + z_rel*Electrons.Get_ElectronTemperature(c+1);
                    break;
                }
            }

            // set ion electron desity and electron temperature
            Neutrals.set_ElectronDensity(i, e_den);
            Neutrals.set_ElectronTemperature(i, e_tmp);
        }
        return Neutrals;
    }

    HypiC::Particles_Object Grid_to_Particles_Ions(HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons){
        // loop over ions
        for(size_t i=0; i<Ions._nParticles; ++i){
            double e_den;
            double e_tmp;
            double E_fld;
            double dz = Electrons.Grid_Step;
            double z_p = Ions.get_Position(i);

            // loop over cells
            for(size_t c=0; c<Electrons._nElectrons; ++c){
                // determine if particle is within cell
                double z_cell = Electrons.Get_CellCenter(c);
                double z_rel = abs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor (s for next cell = z_rel)
                    double s = 1 - z_rel;
                    // calculate partial electron density from this cell and next
                    e_den = s*Electrons.Get_PlasmaDensity(c) + z_rel*Electrons.Get_PlasmaDensity(c+1);
                    // calculate partial electron temperature from this cell and next
                    e_tmp = s*Electrons.Get_ElectronTemperature(c) + z_rel*Electrons.Get_ElectronTemperature(c+1);
                    // calculate partial electric field from this cell and next
                    E_fld = s*Electrons.Get_ElectricField(c) + z_rel*Electrons.Get_ElectricField(c+1);
                    break;
                }
            }

            // set ion electron desity, electron temperature, and electric field
            Ions.set_ElectronDensity(i, e_den);
            Ions.set_ElectronTemperature(i, e_tmp);
            Ions.set_ElectricField(i, E_fld);
        }
        return Ions;
    }
}
