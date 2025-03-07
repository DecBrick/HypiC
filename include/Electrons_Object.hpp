#pragma once
#include <memory>//for ptrs 
#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>
#include "Rate_Table.hpp"
#include "Input.hpp"
#include "math.hpp"

namespace HypiC{
    class Electrons_Object{
        private:

        public:
            Electrons_Object();
            ~Electrons_Object();
            std::vector<double> Cell_Center;
            std::vector<double> Electron_Pressure;
            std::vector<double> Electron_Pressure_Gradient;
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
            std::vector<double> Magnetic_Field_T;
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
            std::vector<double> Source_Energy;
            std::vector<double> Ionization_Rate;
            std::vector<double> Delta_ni; 
            std::vector<bool> Ionize_Flag;

            double Grid_Step = 0.0;
            double Id = 0.0;
            double Delta_ni_sum = 0.0;
            size_t _nElectrons = 0;

            void Add_Electron(double electron_density, double electron_temp,
            double magnetic_field, double energy_density,double electron_velocity, double anom_freq, 
            double Efield, double cell_center);
            void Set_Densities(size_t index, double neutral_density, double plasma_density, double current_density);
            void Set_Velocities(size_t index, double neutral_velocity, double ion_velocity);
            void Update_Delta_ni(size_t index, double dni);
            double Get_CellCenter(size_t index);
            double Get_PlasmaDensity(size_t index);
            double Get_NeutralDensity(size_t index);
            double Get_ElectronTemperature(size_t index);
            double Get_ElectricField(size_t index);
            double Get_Potential(size_t index);
            void Clear_Out_Particles();
            void Update_From_Neutrals(size_t index, double neutral_density, double neutral_velocity, double neutral_energy);
            void Update_From_Ions(size_t index, double plasma_density, double current_density, double ion_velocity);
            void Normalize_Interpolations();

            void Update_Mobility(HypiC::Options_Object Simulation_Parameters, HypiC::Rate_Table_Object Ionization_Rates);

            void Update_Velocity(HypiC::Options_Object Simulation_Parameters);

            void Update_Pressure_Gradient(HypiC::Options_Object Simulation_Parameters);

            void Update_Thermal_Conductivity(HypiC::Options_Object Simulation_Parameters);

            void Compute_Discharge_Current(HypiC::Options_Object Simulation_Parameters);

            void Compute_Electric_Field(HypiC::Options_Object Simulation_Parameters);

            void Solve_Potential(HypiC::Options_Object Simulation_Parameters);

            void Update_Electron_Energy(HypiC::Options_Object Simulation_Parameters, HypiC::Rate_Table_Object Loss_Rates);
    };
}
