#include <iostream> 
#include <cmath>//for exp
#include "Particles_Object.hpp"

namespace HypiC
{

    // Default constructor
    //template <class fp_type>
    //Particles_Object<fp_type>::Particles_Object()
    Particles_Object::Particles_Object()
    {
        std::cout << "Hello from Constructor!\n";
    }
    
    // Destructor
    //template <class fp_type>
    //Particles_Object<fp_type>::~Particles_Object()
    Particles_Object::~Particles_Object()
    {
        std::cout << "Byebye from  Destructor~\n";
    }

    
    //Add particle method
    //template <class fp_type>
    //void Particles_Object<fp_type>::Add_Particle(){
    void Particles_Object::Add_Particle(){
        std::cout << "Hello From Add Particle Method\n";
    }


    //Remove particle method
    //template <class fp_type>
    //void Particles_Object<fp_type>::Remove_Particle(){
    void Particles_Object::Remove_Particle(){
        std::cout << "Hello From Remove Particle Method\n";
    }

    //Update particle method
    //template <class fp_type>
    //void Particles_Object<fp_type>::Update_Particle(size_t index, double dt){
    void Particles_Object::Update_Particle(size_t index, double dt){
        std::cout << "Hello From Update Particle Method\n";
        double rate_coefficient;
        //update the velocity
        this->_Velocities[index]+= this->_ChargetoMassRatio * this->_ElectricField[index] * dt;

        //update the position 
        this->_Positions[index]+= this->_Velocities[index] * dt;

        //update the weights
        //rate_coefficient = Ionization_rate(this->_ElectronTemperature[index]);
        //this->_Weights[index]= this->_Weights[index] + this->_IonizationDirection * exp(-1.0 * this->_ElectronDensity * rate_coefficient * dt);
    }
    
}