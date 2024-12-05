#include <cmath>
#include "HypiCpp.hpp"

namespace HypiC{
    void Particles_to_Grid(HypiC::Particles_Object Neutrals, HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons){
        // loop over cells
        for(size_t c=0; c<Electrons._nElectrons; c++){
            double N_den = 0.0;
            double N_vel = 0.0;
            double P_den = 0.0;
            double I_vel = 0.0;
            double dz = Electrons.Grid_Step;
            double z_cell = 0.5*dz + c*dz;

            // loop over neutrals
            for(size_t i=0; i<Neutrals._nParticles; ++i){
                // determine if particle is within cell
                double z_p = Neutrals.get_Position(i);
                double z_rel = abs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor
                    double s = 1 - z_rel;
                    // weight
                    double w = Neutrals.get_Weight(i);
                    // calculate partial density and velocity
                    N_den += s*w;
                    N_vel += s*w*Neutrals.get_Velocity(i);
                }
            }

            // loop over ions
            for(size_t i=0; i<Ions._nParticles; ++i){
                // determine if particle is within cell
                double z_p = Ions.get_Position(i);
                double z_rel = abs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor
                    double s = 1 - z_rel;
                    // weight
                    double w = Ions.get_Weight(i);
                    // calculate partial density and velocity
                    P_den += s*w;
                    I_vel += s*w*Ions.get_Velocity(i);
                }
            }

            // set cell density and velocity
            Electrons.Set_Densities(c,N_den,P_den);
            Electrons.Set_Velocities(c,N_vel,I_vel);
        }
    }

    void Grid_to_Particles(HypiC::Particles_Object Ions, HypiC::Electrons_Object Electrons){
        // loop over ions
        for(size_t i=0; i<Ions._nParticles; ++i){
            double E_fld = 0.0;
            double dz = Electrons.Grid_Step;
            double z_p = Ions.get_Position(i);

            // loop over cells
            for(size_t c=0; c<Electrons._nElectrons; ++c){
                // determine if particle is within cell
                double z_cell = 0.5*dz + c*dz;
                double z_rel = abs( (z_cell - z_p)/dz);
                if(z_rel < 1){
                    // shape factor
                    double s = 1 - z_rel;
                    // calculate partial electric field
                    E_fld += s*Electrons.Get_ElectricField(c);
                }
            }

            // set ion electric field
            Ions.set_ElectricField(i, E_fld);
        }
    }
}
