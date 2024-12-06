#pragma once
#include <memory>//for ptrs 
#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>
#include "Rate_Table.hpp"

namespace HypiC{
    class Electrons_Object{
        private:

        public:
            Electrons_Object();
            ~Electrons_Object();
            std::vector<double> Cell_Center;
            std::vector<double> Electron_Pressure;
            std::vector<double> Electron_Pressure_Gradient;
            std::vector<double> Electron_Kinetic_Energy;
            std::vector<double> Electron_Thermal_Conductivity;
            std::vector<double> Discharge_Current;
            std::vector<double> EnergyDensity;
            std::vector<double> Neutral_Density_m3;
            std::vector<double> Neutral_Velocity_m_s;
            std::vector<double> Neutral_Temperature_K;
            std::vector<double> Plasma_Density_m3;
            std::vector<double> Ion_Velocity_m_s;
            std::vector<double> Ion_Temperature_eV;
            std::vector<double> Ion_Current_Density;
            std::vector<double> Electron_Velocity_m_s;
            std::vector<double> Electron_Temperature_eV;
            std::vector<double> Magnetic_Field_G;
            std::vector<double> Electric_Field_V_m;
            std::vector<double> Potential;
            std::vector<double> Anomalous_Frequency_Hz;
            std::vector<double> Freq_Elec_Neutral;
            std::vector<double> Freq_Classical;
            std::vector<double> Freq_Electron_Wall_Collision;
            std::vector<double> Freq_Anomalous_Collision;
            std::vector<double> Freq_Total_Electron_Collision;
            std::vector<double> Electron_Mobility;
            std::vector<double> Ion_Z;
            std::vector<double> Ionization_Rate;

            double Grid_Step = 0.0;
            size_t _nElectrons = 0;

            void Add_Electron(double electron_density, double electron_temp,
            double magnetic_field, double energy_density,double electron_velocity, double anom_freq, 
            double Efield, double cell_center);
            void Set_Densities(size_t index, double neutral_density, double plasma_density);
            void Set_Velocities(size_t index, double neutral_velocity, double ion_velocity);
            double Get_ElectricField(size_t index);
    };
}
