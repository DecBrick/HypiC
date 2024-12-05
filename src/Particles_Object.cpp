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
        
    }
    
    // Destructor
    //template <class fp_type>
    //Particles_Object<fp_type>::~Particles_Object()
    Particles_Object::~Particles_Object()
    {
        
    }

    
    //Add particle method
    //template <class fp_type>
    //void Particles_Object<fp_type>::Add_Particle(){
    void Particles_Object::Add_Particle(double init_position, double init_velocity, 
    double init_weight, double init_field, double init_density, double init_temperature){
        //add the values
        this->_Positions.push_back(init_position);
        this->_Velocities.push_back(init_velocity);
        this->_Weights.push_back(init_weight);
        this->_ElectricField.push_back(init_field);
        this->_ElectronDensity.push_back(init_density);
        this->_ElectronTemperature.push_back(init_temperature);
        //increase n particles count
        this->_nParticles+=1;
    }


    //Remove particle method
    //template <class fp_type>
    //void Particles_Object<fp_type>::Remove_Particle(){
    void Particles_Object::Remove_Particle(size_t index){
        assert(this->_nParticles > 0);//check that we've initialized before
        //remove values from the vectors
        this->_Positions.erase(this->_Positions.begin() + index);
        this->_Velocities.erase(this->_Velocities.begin() + index);
        this->_Weights.erase(this->_Weights.begin() + index);
        this->_ElectricField.erase(this->_ElectricField.begin() + index);
        this->_ElectronDensity.erase(this->_ElectronDensity.begin() + index);
        this->_ElectronTemperature.erase(this->_ElectronTemperature.begin() + index);

        //update the number of particles
        this->_nParticles-=1;
    }

    //Update particle method
    void Particles_Object::Update_Particle(size_t index, double dt, HypiC::Rate_Table_Object Single_Ionization_Table){
        double rate_coefficient;
        //update the velocity
        this->_Velocities[index]+= this->_ChargetoMassRatio * this->_ElectricField[index] * dt;

        //update the position 
        this->_Positions[index]+= this->_Velocities[index] * dt;

        //update the weights
        rate_coefficient = Single_Ionization_Table.interpolate(this->_ElectronTemperature[index]);
        this->_Weights[index] *= exp(this->_IonizationDirection * this->_ElectronDensity[index] * rate_coefficient * dt);
    }

    void Particles_Object::Velocity_Backstep(double dt){
        //based on the current electric field, set v back by dt/2 to initialize leapfrog
        for (size_t i=0; i<this->_nParticles; ++i){
            this->_Velocities[i] -= (dt/2) * this->_ChargetoMassRatio * this->_ElectricField[i];
        }
    }
    
    void Particles_Object::set_Position(size_t index, double value){
        this->_Positions[index] = value;
    }
    void Particles_Object::set_Velocity(size_t index, double value){
        this->_Velocities[index] = value;
    }

    void Particles_Object::set_ElectricField(size_t index, double value){
        this->_ElectricField[index] = value;
    }

    //accessor methods
    double Particles_Object::get_Position(size_t index){
        return this->_Positions[index];
    }
    double Particles_Object::get_Velocity(size_t index){
        return this->_Velocities[index];
    }
    double Particles_Object::get_Weight(size_t index){
        return this->_Weights[index];
    }
    double Particles_Object::get_ElectricField(size_t index){
        return this->_ElectricField[index];
    }
    double Particles_Object::get_ElectronDensity(size_t index){
        return this->_ElectronDensity[index];
    }
    double Particles_Object::get_ElectronTemperature(size_t index){
        return this->_ElectronTemperature[index];
    }
}