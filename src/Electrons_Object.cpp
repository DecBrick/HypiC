#include <iostream> 
#include <cmath>//for exp
#include "Electrons_Object.hpp"

namespace HypiC
{

    // Default constructor
    //template <class fp_type>
    //Particles_Object<fp_type>::Particles_Object()
    Electrons_Object::Electrons_Object()
    {
        
    }
    
    // Destructor
    //template <class fp_type>
    //Particles_Object<fp_type>::~Particles_Object()
    Electrons_Object::~Electrons_Object()
    {      
    }

    void Electrons_Object::Add_Electron(double ne, double ChargeDensity){
        //add the values
        this->_ElecDensity.push_back(ne);
        this->_ChargeDensity.push_back(ChargeDensity);
        //increase n particles count
        this->_nElectrons+=1;
    }

}