#include <cmath>
#include <iostream>
#include "HypiCpp.hpp"

namespace HypiC{
    HypiC::Electrons_Object Particles_to_Grid(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons){
        // loop over cells
        for(size_t c=0; c<Electrons._nElectrons; ++c){
            double N_den = 0.0;
            double N_flux = 0.0;
            double P_den = 0.0;
            double I_vel = 0.0;
            double J_den = 0.0; // current density
            double dz = Electrons.Grid_Step;
            double z_cell = Electrons.Cell_Center[c];

            // loop over neutrals
            for(size_t i=0; i<Neutrals._nParticles; ++i){
                // determine if particle is within cell
                double z_p = Neutrals.get_Position(i);
                double z_rel = abs( (z_cell - z_p)/(dz));
                if(z_rel < 1){
                    // shape factor
                    double s = 1 - z_rel;
                    // weight
                    double w = Neutrals.get_Weight(i);
                    // calculate partial number density
                    N_den += s*w /dz;
                    // calculaye partial velocity
                    N_flux += s*(w/dz)*Neutrals.get_Velocity(i);
                }
            }

            // loop over ions
            for(size_t i=0; i<Ions._nParticles; ++i){
                // determine if particle is within cell
                double z_p = Ions.get_Position(i);
                double z_rel = abs( (z_cell - z_p)/(dz));
                if(z_rel < 1){
                    // shape factor
                    double s = 1 - z_rel;
                    // weight
                    double w = Ions.get_Weight(i);
                    // calculate partial number density
                    P_den += s*w / dz;
                    // calculate partial current density
                    J_den += s*(w/dz)*Ions.get_Velocity(i);
                    /*if (c == 5) { 
                        std::cout << "---------------------\n";
                        std::cout << i <<"\n";
                        std::cout << "---------------------\n";
                        std::cout << Ions.get_Position(i) <<"\n";
                        std::cout << Ions.get_Weight(i) <<"\n";
                        std::cout << Ions.get_Velocity(i) <<"\n";
                    }*/
                }
            }

            // set cell densities and velocities
            Electrons.Set_Densities(c,N_den,P_den,1.602176634e-19 * J_den);
            Electrons.Set_Velocities(c,N_flux / N_den ,J_den / P_den);
        }

        return Electrons;

    }

    HypiC::Particles_Object Grid_to_Particles_Neutrals(HypiC::Particles_Object Neutrals, HypiC::Electrons_Object Electrons){
        // loop over neutrals
        for(size_t i=0; i<Neutrals._nParticles; ++i){
            double e_den = 0.0;
            double e_tmp = 0.0;
            double dz = Electrons.Grid_Step;
            double z_p = Neutrals.get_Position(i);

            // loop over cells
            for(size_t c=0; c<Electrons._nElectrons; ++c){
                // determine if particle is within cell
                double z_cell = Electrons.Get_CellCenter(c);
                double z_rel = abs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor
                    double s = 1 - z_rel;
                    // calculate partial electron density
                    e_den += s*Electrons.Get_PlasmaDensity(c);
                    // calculate partial electron temperature
                    e_tmp += s*Electrons.Get_ElectronTemperature(c);
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
            double e_den = 0.0;
            double e_tmp = 0.0;
            double E_fld = 0.0;
            double dz = Electrons.Grid_Step;
            double z_p = Ions.get_Position(i);

            // loop over cells
            for(size_t c=0; c<Electrons._nElectrons; ++c){
                // determine if particle is within cell
                double z_cell = Electrons.Get_CellCenter(c);
                double z_rel = abs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor
                    double s = 1 - z_rel;
                    // calculate partial electron density
                    e_den += s*Electrons.Get_PlasmaDensity(c);
                    // calculate partial electron temperature
                    e_tmp += s*Electrons.Get_ElectronTemperature(c);
                    // calculate partial electric field
                    E_fld += s*Electrons.Get_ElectricField(c);
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
