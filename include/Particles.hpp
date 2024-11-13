#pragma once
#include <HypiCpp.hpp>
#include <memory>


namespace HypiC{

    template <class fp_type>
    class Particles_Object{
        private: 
            size_t _nParticles = 0;
            std::unique_ptr< size_t [] > _EmptyIdx = nullptr;
        public: 
            std::unique_ptr< fp_type [] > _Positions = nullptr;
            std::unique_ptr< fp_type [] > _Velocities = nullptr;
            std::unique_ptr< fp_type [] > _Weights = nullptr;
            std::unique_ptr< fp_type [] > _ElectricField = nullptr;
            std::unique_ptr< fp_type [] > _ElectronDensity = nullptr;
            std::unique_ptr< fp_type [] > _ElectronTemperature = nullptr;
            fp_type _ChargetoMassRatio = 0;
            fp_type _IonizationDirection = 0;//flag, should be -1 for neutrals, +1 for ions

            Particles_Object();
            virtual ~Particles_Object();
            void Add_Particle();
            void Remove_Particle();
            void Update_Particle(size_t index, double dt);
    };
}