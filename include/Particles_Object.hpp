#pragma once
#include <memory>//for ptrs 
#include <cassert>
#include <cstddef>
#include <map>
#include <utility>



namespace HypiC{

    //template <class fp_type>
    class Particles_Object{
        private: 
            size_t _nParticles = 0;
            std::unique_ptr< size_t [] > _EmptyIdx = nullptr;
        public: 
            std::unique_ptr< double [] > _Positions = nullptr;
            std::unique_ptr< double [] > _Velocities = nullptr;
            std::unique_ptr< double [] > _Weights = nullptr;
            std::unique_ptr< double [] > _ElectricField = nullptr;
            std::unique_ptr< double [] > _ElectronDensity = nullptr;
            std::unique_ptr< double [] > _ElectronTemperature = nullptr;
            double _ChargetoMassRatio = 0;
            int _IonizationDirection = 0;//flag, should be -1 for neutrals, +1 for ions

            Particles_Object();
            ~Particles_Object();
            void Add_Particle();
            void Remove_Particle();
            void Update_Particle(size_t index, double dt);
    };
}