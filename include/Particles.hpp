#pragma once
#include <HypiCpp.hpp>

namespace HypiC{

    template <class fp_type>
    class Particles_Object{
        private: 
            size_t _nParticles = 0;
            std::unique_ptr< size_t [] > _EmptyIdx = nullptr;
        public: 
            std::unique_ptr< fp_type [] > _Positions = nullptr;
            std::unique_ptr< fp_type [] > _Velocities = nullptr;
            std::unique_ptr< fp_type [] > _ElectricField = nullptr;
            std::unique_ptr< fp_type [] > _ElectronDensity = nullptr;
            std::unique_ptr< fp_type [] > _ElectronTemperature = nullptr;

            Particles_Object();
            virtual ~Particles_Object();
            void Add_Particle();
            void Remove_Particle();
            void Update_Particle();
    };
}