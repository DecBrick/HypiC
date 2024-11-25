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

    void Electrons_Object::Add_Electron(double electron_density, double electron_temp,
    double charge_density, double magnetic_field){
        //add the values
        this->Electric_Field_V_m.push_back(electron_density);
        this->Electron_Temperature_eV.push_back(electron_temp);
        this->EnergyDensity.push_back(charge_density);
        this->Magnetic_Field_G.push_back(magnetic_field);
        //increase n particles count
        this->_nElectrons+=1;
    }

}