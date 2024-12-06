#pragma once
#include <memory>//for ptrs 
#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>
#include "Rate_Table.hpp"



namespace HypiC{

    //template <class fp_type>
    class Particles_Object{
        private: 
            
        public: 
            size_t _nParticles = 0;
            std::vector<double> _Positions;
            std::vector<double> _Velocities;
            std::vector<double> _Weights;
            std::vector<double> _ElectricField;
            std::vector<double> _ElectronDensity;
            std::vector<double> _ElectronTemperature;
            double _ChargetoMassRatio = 0;
            int _IonizationDirection = 0;//flag, should be -1 for neutrals, +1 for ions
            //construction/destruction methods
            Particles_Object();
            ~Particles_Object();
            //manipulating particle methods
            void Add_Particle(double init_position, double init_velocity, 
            double init_weight, double init_field, double init_density, double init_temperature);
            void Remove_Particle(size_t index);
            void Update_Particle(size_t index, double dt, HypiC::Rate_Table_Object);
            void set_Position(size_t index, double value);
            void set_Velocity(size_t index, double value);
            void set_ElectricField(size_t index, double value);
            void Velocity_Backstep(double dt);
            //accessor methods
            double get_Position(size_t index);
            double get_Velocity(size_t index);
            double get_Weight(size_t index);
            double get_ElectricField(size_t index);
            double get_ElectronDensity(size_t index);
            double get_ElectronTemperature(size_t index);

    };
}