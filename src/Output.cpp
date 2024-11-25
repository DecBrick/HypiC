#include "HypiCpp.hpp"
#include <fstream>
#include <iomanip>

namespace HypiC{
    void Time_Sum_Object::Initialize_Time_Sum(size_t nCells){
        //set time equal to 0
        this->time = 0;
        //set equal to the size of the grid 
        for (size_t i=0; i<nCells; ++i){
            this->Neutral_Density_m3.push_back(0.0);
            this->Neutral_Velocity_m_s.push_back(0.0);
            this->Neutral_Temperature_K.push_back(0.0);
            this->Plasma_Density_m3.push_back(0.0);
            this->Ion_Velocity_m_s.push_back(0.0);
            this->Ion_Temperature_eV.push_back(0.0);
            this->Electron_Velocity_m_s.push_back(0.0);
            this->Electron_Temperature_eV.push_back(0.0);
            this->Magnetic_Field_G.push_back(0.0);
            this->Electric_Field_V_m.push_back(0.0);
            this->Anomalous_Frequency_Hz.push_back(0.0);
        }
    };
    //might have to add some averaging here as the interpolation step might not consider temperature
    /*void Time_Sum_Object::Time_Sum(HypiC::Electron_Object Quantities){
        //increment time
        time += Quantities.dt
        //sum quantities
        for (size_t i=0; i<Inputs.nCells; ++i){
            this->Neutral_Density_m3[i] += Electron_Object.nn[i];
            this->Neutral_Velocity_m_s[i] += Electron_Object.nv[i];
            this->Neutral_Temperature_K[i] += Electron_Object.Tn[i];
            this->Plasma_Density_m3[i] += Electron_Object.ne[i];
            this->Ion_Velocity_m_s[i] += Electron_Object.vz[i];
            this->Ion_Temperature_eV[i] += Electron_Object.Ti[i];
            this->Electron_Velocity_m_s[i] += Electron_Object.ev[i];
            this->Electron_Temperature_eV[i] += Electron_Object.Te[i];
            this->Magnetic_Field_G[i] += Electron_Object.B[i];
            this->Electric_Field_V_m[i] += Electron_Object.E[i];
            this->Anomalous_Frequency_Hz[i] += Electron_Object.nu_an[i];
        }
    };*/

    void Time_Sum_Object::Write_Output(std::string Filename, size_t nCells){
        int digits = 3;
        //write to csv
        //open the file 
        std::ofstream f;
        f.open(Filename);
        //write the headers
        f << "nn_m3, nv_m_s, Tn_K, ne_m3, vi_m_s, Ti_eV, ve_m_s, Te_eV, B_G, E_V_m, nu_an_Hz\n";
        //write the quantities of interest
        for (size_t i=0; i<nCells; ++i){
            f << std::setprecision(digits) << this->Neutral_Density_m3[i] / this->time << ",";
            f << std::setprecision(digits) << this->Neutral_Velocity_m_s[i] / this->time << ",";
            f << std::setprecision(digits) << this->Neutral_Temperature_K[i] / this->time << ",";
            f << std::setprecision(digits) << this->Plasma_Density_m3[i] / this->time << ",";
            f << std::setprecision(digits) << this->Ion_Velocity_m_s[i] / this->time << ",";
            f << std::setprecision(digits) << this->Ion_Temperature_eV[i] / this->time << ",";
            f << std::setprecision(digits) << this->Electron_Velocity_m_s[i] / this->time << ",";
            f << std::setprecision(digits) << this->Electron_Temperature_eV[i] / this->time << ",";
            f << std::setprecision(digits) << this->Magnetic_Field_G[i] / this->time << ",";
            f << std::setprecision(digits) << this->Electric_Field_V_m[i] / this->time << ",";
            f << std::setprecision(digits) << this->Anomalous_Frequency_Hz[i] / this->time << ",";
            f << "\n";
        }

        //close the file
        f.close();
    };
}