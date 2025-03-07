#include "HypiCpp.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace HypiC{
    void Time_Sum_Object::Initialize_Time_Sum(size_t nCells, HypiC::Electrons_Object Grid){
        //set time and discharge current equal to 0
        this->time = 0.0;
        this->Discharge_Current =0.0;
        this->Id_Last = 0.0;
        //set equal to the size of the grid 
        for (size_t i=0; i<nCells; ++i){
            this->z_m.push_back(Grid.Cell_Center[i]);
            this->Neutral_Density_m3.push_back(0.0);
            this->Neutral_Velocity_m_s.push_back(0.0);
            this->Plasma_Density_m3.push_back(0.0);
            this->Ion_Velocity_m_s.push_back(0.0);
            this->Electron_Velocity_m_s.push_back(0.0);
            this->Electron_Temperature_eV.push_back(0.0);
            this->Magnetic_Field_G.push_back(Grid.Magnetic_Field_T[i] * 10000.0);
            this->Electric_Field_V_m.push_back(0.0);
            this->Anomalous_Frequency_Hz.push_back(0.0);
            this->Ionization_Rate_m3_s.push_back(0.0);
            this->Potential_V.push_back(0.0);
        }
    };
    //might have to add some averaging here as the interpolation step might not consider temperature
    void Time_Sum_Object::Time_Sum(HypiC::Electrons_Object Electrons, HypiC::Options_Object Simulation_Parameters){
        //increment time
        this->time += Simulation_Parameters.dt;
        
        //pull timestep
        double dt = Simulation_Parameters.dt;

        //sum quantities
        this->Discharge_Current += dt * Electrons.Id;
        this->Id_Last = Electrons.Id;
        //#pragma omp parallel for
        for (size_t i=0; i<Simulation_Parameters.nCells; ++i){
            this->Neutral_Density_m3[i] += dt*Electrons.Neutral_Density_m3[i];
            this->Neutral_Velocity_m_s[i] += dt*Electrons.Neutral_Velocity_m_s[i];
            this->Plasma_Density_m3[i] += dt*Electrons.Plasma_Density_m3[i];
            this->Ion_Velocity_m_s[i] += dt*Electrons.Ion_Velocity_m_s[i];
            this->Electron_Velocity_m_s[i] += dt*Electrons.Electron_Velocity_m_s[i];
            this->Electron_Temperature_eV[i] += dt*Electrons.Electron_Temperature_eV[i];
            this->Electric_Field_V_m[i] += dt*Electrons.Electric_Field_V_m[i];
            this->Anomalous_Frequency_Hz[i] += dt*Electrons.Anomalous_Frequency_Hz[i];
            this->Ionization_Rate_m3_s[i] += dt*Electrons.Ionization_Rate[i];
            this->Potential_V[i] += dt*Electrons.Potential[i];
        }
    };

    void Time_Sum_Object::Write_Output(std::string Filename, size_t nCells){
        int digits = 3;
        //for unit testing
        if (this->time == 0){
            this->time = 1;
        }
        //write to csv
        //open the file 
        std::ofstream f;
        f.open(Filename);
        //write the headers
        f << "z_m, nn_m3, nv_m_s, ne_m3, vi_m_s, ve_m_s, Te_eV, B_G, E_V_m, nu_an_Hz, k_iz_m-3_s-1, phi_V\n";
        //write the quantities of interest
        for (size_t i=0; i<nCells; ++i){
            f << std::setprecision(digits) << this->z_m[i] << ",";
            f << std::setprecision(digits) << this->Neutral_Density_m3[i] / this->time << ",";
            f << std::setprecision(digits) << this->Neutral_Velocity_m_s[i] / this->time << ",";
            f << std::setprecision(digits) << this->Plasma_Density_m3[i] / this->time << ",";
            f << std::setprecision(digits) << this->Ion_Velocity_m_s[i] / this->time << ",";
            f << std::setprecision(digits) << this->Electron_Velocity_m_s[i] / this->time << ",";
            f << std::setprecision(digits) << this->Electron_Temperature_eV[i] / this->time << ",";
            f << std::setprecision(digits) << this->Magnetic_Field_G[i] << ",";
            f << std::setprecision(digits) << this->Electric_Field_V_m[i] / this->time << ",";
            f << std::setprecision(digits) << this->Anomalous_Frequency_Hz[i] / this->time << ",";
            f << std::setprecision(digits) << this->Ionization_Rate_m3_s[i] / this->time << ",";
            f << std::setprecision(digits) << this->Potential_V[i] / this->time << ",";
            f << "\n";
        }

        //close the file
        f.close();
        std::cout << "Simulation Time is: " << time << " Discharge Current is: " << this->Id_Last << " Average Discharge Current is: "<< this->Discharge_Current/this->time << "\n";
        std::cout << "-----------------------------------\n";
    };
}