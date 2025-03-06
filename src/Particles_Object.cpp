#include <iostream> 
#include <cmath>//for exp
#include "Particles_Object.hpp"
#include "HypiCpp.hpp"

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
    double init_weight, size_t init_cell, int init_cell2, double init_field){
        //add the values
        this->_Positions.push_back(init_position);
        this->_Velocities.push_back(init_velocity);
        this->_Weights.push_back(init_weight);
        this->_CellIndex.push_back(init_cell);
        this->_Cell2Index.push_back(init_cell2);
        this->_ElectricField.push_back(init_field);
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
        this->_CellIndex.erase(this->_CellIndex.begin() + index);
        this->_Cell2Index.erase(this->_Cell2Index.begin() + index);
        this->_ElectricField.erase(this->_ElectricField.begin() + index);

        //update the number of particles
        this->_nParticles-=1;

        //try setting properites to 0
        /*this->_Positions[index] = 0.0;
        this->_Velocities[index] = 0.0;
        this->_Weights[index] = 0.0;*/
    }

    //Update particle method
    void Particles_Object::Push_Particle(size_t index, double dt, HypiC::Electrons_Object Electrons){
        double rate_coefficient;
        //update the velocity
        this->_Velocities[index]+= this->_ChargetoMassRatio * this->_ElectricField[index] * dt;

        //update the position 
        this->_Positions[index]+= this->_Velocities[index] * dt;

        //update the cell
        //if the distance is more than half a cell width away from the cell center
        if (fabs(this->_Positions[index] - Electrons.Cell_Center[this->_CellIndex[index]]) > (Electrons.Grid_Step / 2.0)){
            //adjust the cell number, the velocity comparision is to account for the direction the particle has traveled
            //the static cast accounts for the number of cells it has traveled 
            this->_CellIndex[index] += ((this->_Velocities[index] > 0.0) - (this->_Velocities[index] < 0.0)) * 
            static_cast<size_t>((fabs(this->_Positions[index] - Electrons.Cell_Center[this->_CellIndex[index]]) / Electrons.Grid_Step)+0.5);
        }
        //update cell2 index (either +-1 dependening on relative position), accounts for new cell 
        if ((this->_Positions[index] - Electrons.Cell_Center[this->_CellIndex[index]]) * this->_Cell2Index[index] < 0){
            this->_Cell2Index[index] *=-1;//flip the sign
        }
    }

    void Particles_Object::Velocity_Backstep(double dt){
        //based on the current electric field, set v back by dt/2 to initialize leapfrog
        #pragma omp parallel for
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
    void Particles_Object::set_Weight(size_t index, double value){
        this->_Weights[index] = value;
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
    size_t Particles_Object::get_Cell(size_t index){
        return this->_CellIndex[index];
    }
    int Particles_Object::get_Cell2(size_t index){
        return this->_Cell2Index[index];
    }
    double Particles_Object::get_ElectricField(size_t index){
        return this->_ElectricField[index];
    }
}