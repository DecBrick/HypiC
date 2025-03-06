#pragma once
#include <memory>//for ptrs 
#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>
#include "Rate_Table.hpp"
#include "Electrons_Object.hpp"



namespace HypiC{

    //template <class fp_type>
    class Particles_Object{
        private: 
            
        public: 
            size_t _nParticles = 0;
            std::vector<double> _Positions;
            std::vector<double> _Velocities;
            std::vector<double> _Weights;
            std::vector<size_t> _CellIndex;
            std::vector<int> _Cell2Index;
            std::vector<double> _ElectricField;
            std::vector<bool> _Ionization_Flag;
            double _ChargetoMassRatio = 0;
            //construction/destruction methods
            Particles_Object();
            ~Particles_Object();
            //manipulating particle methods
            void Add_Particle(double init_position, double init_velocity, 
            double init_weight, size_t init_cell, int init_cell2, double init_field);
            void Remove_Particle(size_t index);
            void Push_Particle(size_t index, double dt, HypiC::Electrons_Object);
            void set_Position(size_t index, double value);
            void set_Velocity(size_t index, double value);
            void set_Weight(size_t index, double value);
            void set_ElectricField(size_t index, double value);
            void Velocity_Backstep(double dt);
            //accessor methods
            double get_Position(size_t index);
            double get_Velocity(size_t index);
            double get_Weight(size_t index);
            size_t get_Cell(size_t index);
            int get_Cell2(size_t index);
            double get_ElectricField(size_t index);

    };
}