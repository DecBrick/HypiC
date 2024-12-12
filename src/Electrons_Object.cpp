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
        this->Electron_Kinetic_Energy.push_back(1.5 * electron_density * electron_temp);
        this->Electron_Pressure.push_back(1.602176634e-19 * electron_density * electron_temp);
        this->Magnetic_Field_T.push_back(magnetic_field);
        this->EnergyDensity.push_back(energy_density);
        this->Electron_Velocity_m_s.push_back(electron_velocity);
        this->Anomalous_Frequency_Hz.push_back(anom_freq);
        this->Electric_Field_V_m.push_back(Efield);
        
        //set other inital properties to 0 to correctly initalize vector size
        this->Neutral_Density_m3.push_back(0.0);
        this->Ion_Current_Density.push_back(0.0);
        this->Freq_Elec_Neutral.push_back(0.0);
        this->Freq_Classical.push_back(0.0);
        this->Freq_Anomalous_Collision.push_back(0.0);
        this->Freq_Electron_Wall_Collision.push_back(0.0);
        this->Freq_Total_Electron_Collision.push_back(0.0);
        this->Ionization_Rate.push_back(0.0);
        this->Electron_Mobility.push_back(0.0);
        this->Electron_Pressure_Gradient.push_back(0.0);
        this->Potential.push_back(0.0);
        this->Electron_Thermal_Conductivity.push_back(0.0);
        this->Neutral_Velocity_m_s.push_back(0.0);
        this->Neutral_Temperature_K.push_back(0.0);
        this->Ion_Velocity_m_s.push_back(0.0);
        this->Ion_Temperature_eV.push_back(0.0);


        //increase n particles count
        this->_nElectrons+=1;
    }

    void Electrons_Object::Set_Densities(size_t index, double neutral_density, double plasma_density, double current_density){
        this->Neutral_Density_m3[index] = neutral_density;
        this->Plasma_Density_m3[index] = plasma_density;
        this->Ion_Current_Density[index] = current_density;
    }

    void Electrons_Object::Set_Velocities(size_t index, double neutral_velocity, double ion_velocity){
        this->Neutral_Velocity_m_s[index] = neutral_velocity;
        this->Ion_Velocity_m_s[index] = ion_velocity;
    }

    double Electrons_Object::Get_CellCenter(size_t index){
        return this->Cell_Center[index];
    }
    
    double Electrons_Object::Get_PlasmaDensity(size_t index){
        return this->Plasma_Density_m3[index];
    }

    double Electrons_Object::Get_ElectronTemperature(size_t index){
        return this->Electron_Temperature_eV[index];
    }

    double Electrons_Object::Get_ElectricField(size_t index){
        return this->Electric_Field_V_m[index];
    }


    void Electrons_Object::Update_Mobility(HypiC::Options_Object Simulation_Parameters, HypiC::Rate_Table_Object Ionization_Rates){
        double Elec_Cycl_Freq;
        double Beta;
        double Omega;
        
        for(size_t c=0; c<this->_nElectrons; ++c){
            this->Electron_Temperature_eV[c] = (2.0/3.0) * this->EnergyDensity[c] /this->Plasma_Density_m3[c];
            //P = n * e * T
            this->Electron_Pressure[c] =  (2.0/3.0) * 1.602176634e-19 * this->EnergyDensity[c];

        
            //This seems to be super complicated?? But Landmark seems to be constant 2.5e-13 rate?
            this->Freq_Elec_Neutral[c] = 2.5e-13 * this->Neutral_Density_m3[c];
            
            //Electrons.Freq_Classical[c] = Electrons.Freq_Elec_Neutral[c] + Electrons.Freq_Elec_Ion[c];
            this->Freq_Classical[c] = this->Freq_Elec_Neutral[c];

            //Seems like Radial Loss Freq is a constant 1e7 for Landmark?
            if (this->Cell_Center[c] <= Simulation_Parameters.Channel_Length_m) {
                this->Freq_Electron_Wall_Collision[c] = 1e7 ;
            }
            //update the ionization rate
            this->Ionization_Rate[c] = Ionization_Rates.interpolate(this->Electron_Temperature_eV[c]);
            //update the anomalous frequency 
            Elec_Cycl_Freq = 1.602176634e-19 * this->Magnetic_Field_T[c] / 9.10938356e-31;

            //Might need to change this to 0.1, 1?
            if (this->Cell_Center[c] <= Simulation_Parameters.Channel_Length_m){
               Beta = 0.1 / 16.0; 
            }else{
               Beta = 1.0 / 16.0;
            }

            this->Freq_Anomalous_Collision[c] = Elec_Cycl_Freq * Beta;


            //update total collision frequency
            this->Freq_Total_Electron_Collision[c] = this->Freq_Electron_Wall_Collision[c] + this->Freq_Anomalous_Collision[c] + this->Freq_Classical[c];

            Omega =  1.602176634e-19 * this->Magnetic_Field_T[c] / (9.10938356e-31 *this->Freq_Total_Electron_Collision[c]);
            this->Electron_Mobility[c] = 1.602176634e-19 / (9.10938356e-31 * this->Freq_Total_Electron_Collision[c] * (1+pow(Omega,2.0)));
        }
    }

    void Electrons_Object::Update_Pressure_Gradient(HypiC::Options_Object Simulation_Parameters){
        
        this->Electron_Pressure_Gradient[0] = HypiC::Forward_Difference(this->Electron_Pressure[0],this->Electron_Pressure[1],this->Cell_Center[0],this->Cell_Center[1]);
        for(size_t i=1; i<Simulation_Parameters.nCells-1; ++i){
            this->Electron_Pressure_Gradient[i] = HypiC::Central_Difference(this->Electron_Pressure[i-1],this->Electron_Pressure[i+1],this->Cell_Center[i-1],this->Cell_Center[i+1]);
        }
        size_t end = this->Electron_Pressure_Gradient.size() - 1;
        this->Electron_Pressure_Gradient[this->_nElectrons-1] = HypiC::Backward_Difference(this->Electron_Pressure[this->_nElectrons-2],this->Electron_Pressure[this->_nElectrons-1],this->Cell_Center[this->_nElectrons-2],this->Cell_Center[this->_nElectrons-1]);
    }
}