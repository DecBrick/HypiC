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
            std::vector<double> EnergyDensity;
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

            size_t _nElectrons = 0;

        void Add_Electron(double electron_density, double electron_temp,
        double magnetic_field, double energy_density,double electron_velocity, double anom_freq, 
        double Efield);
    };
}
