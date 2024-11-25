#pragma once
#include <vector>
#include <string>

namespace HypiC{
    
    class Time_Sum_Object{
        private:
            double time;
            std::vector<double> Neutral_Density_m3;
            std::vector<double> Neutral_Velocity_m_s;
            std::vector<double> Neutral_Temperature_K;
            std::vector<double> Plasma_Density_m3;
            std::vector<double> Ion_Velocity_m_s;
            std::vector<double> Ion_Temperature_eV;
            std::vector<double> Electron_Velocity_m_s;
            std::vector<double> Electron_Temperature_eV;
            std::vector<double> Magnetic_Field_G;
            std::vector<double> Electric_Field_V_m;
            std::vector<double> Anomalous_Frequency_Hz;
        public:
            void Initialize_Time_Sum(size_t nCells);
            //void Time_Sum(HypiC::Electron_Object Quantities);
            void Write_Output(std::string Filename, size_t nCells);

    };

}