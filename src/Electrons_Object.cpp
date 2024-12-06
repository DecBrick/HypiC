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
    double magnetic_field, double energy_density,double electron_velocity, double anom_freq, 
    double Efield, double cell_center){
        //add the values
        this->Cell_Center.push_back(cell_center);
        this->Plasma_Density_m3.push_back(electron_density);
        this->Electron_Temperature_eV.push_back(electron_temp);
        this->Magnetic_Field_G.push_back(magnetic_field);
        this->EnergyDensity.push_back(energy_density);
        this->Electron_Velocity_m_s.push_back(electron_velocity);
        this->Anomalous_Frequency_Hz.push_back(anom_freq);
        this->Electric_Field_V_m.push_back(Efield);
        
        //increase n particles count
        this->_nElectrons+=1;
    }

    void Electrons_Object::Set_Densities(size_t index, double neutral_density, double plasma_density){
        this->Neutral_Density_m3[index] = neutral_density;
        this-> Plasma_Density_m3[index] = plasma_density;
    }

    void Electrons_Object::Set_Velocities(size_t index, double neutral_velocity, double ion_velocity){
        this->Neutral_Velocity_m_s[index];
        this->Ion_Velocity_m_s[index];
    }

    double Electrons_Object::Get_ElectricField(size_t index){
        return this->Electric_Field_V_m[index];
    }

}