#include <iostream> 
#include <cmath>//for exp
#include "Electrons_Object.hpp"
#include <stdexcept>
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

    void Electrons_Object::Clear_Out_Particles(size_t nCells){
        this->Neutral_Density_m3.clear();
        this->Neutral_Density_m3.resize(nCells, 0.0);
        this->Plasma_Density_m3.clear();
        this->Plasma_Density_m3.resize(nCells, 0.0);
        this->Ion_Current_Density.clear();
        this->Ion_Current_Density.resize(nCells, 0.0);
        this->Neutral_Velocity_m_s.clear();
        this->Neutral_Velocity_m_s.resize(nCells, 0.0);
        this->Ion_Velocity_m_s.clear();
        this->Ion_Velocity_m_s.resize(nCells, 0.0);
    }

    void Electrons_Object::Update_From_Neutrals(size_t index, double neutral_density, double neutral_velocity){
        this->Neutral_Density_m3[index] += neutral_density;
        this->Neutral_Velocity_m_s[index] += neutral_velocity;
    }

    void Electrons_Object::Update_From_Ions(size_t index, double plasma_density, double current_density, double ion_velocity){
        this->Plasma_Density_m3[index] += plasma_density;
        this->Ion_Current_Density[index] += current_density;
        this->Ion_Velocity_m_s[index] += ion_velocity;
    }

    void Electrons_Object::Normalize_Velocities(){
        //for each cell
        for (size_t i = 0; i < this->_nElectrons; ++i){
            if (this->Plasma_Density_m3[i] <= 0){
                //std::cout << "Say Something\n";
            }
            //divide neutral flux sum by total density
            this->Neutral_Velocity_m_s[i] /= this->Neutral_Density_m3[i];
            //divide ion flux sum by total density
            this->Ion_Velocity_m_s[i] /= this->Plasma_Density_m3[i];
            //apply charge factor to current density
            this->Ion_Current_Density[i] *= 1.602176634e-19;
        }
    }

    void Electrons_Object::Update_Mobility(HypiC::Options_Object Simulation_Parameters, HypiC::Rate_Table_Object Ionization_Rates){
        double Elec_Cycl_Freq;
        double Beta;
        double Omega;
        
        for(size_t c=0; c<this->_nElectrons; ++c){
            this->Electron_Temperature_eV[c] = (2.0/3.0) * this->EnergyDensity[c] /this->Plasma_Density_m3[c];
            if(this->Electron_Temperature_eV[c] > 200){
                std::cout << c << "\n";
                std::cout << this->Electron_Temperature_eV[c] << "\n";
                std::cout << this->Plasma_Density_m3[c] << "\n";
                //throw std::invalid_argument("Temperature too high");
            }
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
            this->Ionization_Rate[c] = Ionization_Rates.interpolate(1.5 * this->Electron_Temperature_eV[c]);
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

    void Electrons_Object::Update_Velocity(HypiC::Options_Object Simulation_Parameters){
        for(size_t c1=0; c1<Simulation_Parameters.nCells; ++c1){
            this->Electron_Velocity_m_s[c1] = (this->Ion_Current_Density[c1] - this->Id/Simulation_Parameters.Channel_Area_m2/1.602176634e-19 / this->EnergyDensity[c1]);
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

    void Electrons_Object::Update_Thermal_Conductivity(HypiC::Options_Object Simulation_Parameters){
        for(size_t i=0; i<Simulation_Parameters.nCells; ++i){
            this->Electron_Thermal_Conductivity[i] = (10.0 / (9.0)) * this->Electron_Mobility[i] * this->EnergyDensity[i];
        }
    }

    void Electrons_Object::Compute_Discharge_Current(HypiC::Options_Object Simulation_Parameters){
        double int1 = 0;
        double int2 = 0;
        double Dz = Simulation_Parameters.Domain_Length_m / Simulation_Parameters.nCells;
        double dz;
        double int1_1;
        double int1_2;
        double int2_1;
        double int2_2;
        double enemu1;
        double enemu2;

        for (size_t c=0; c<Simulation_Parameters.nCells-1; ++c){
            //dz = Dz * c;
            enemu1 = 1.602176634e-19 * this->Plasma_Density_m3[c] * this->Electron_Mobility[c];
            enemu2 = 1.602176634e-19 * this->Plasma_Density_m3[c+1] * this->Electron_Mobility[c+1];
            int1_1 = (this->Ion_Current_Density[c]/enemu1 + this->Electron_Pressure_Gradient[c]/this->Plasma_Density_m3[c]/1.602176634e-19 );
            int1_2 = (this->Ion_Current_Density[c+1]/enemu2 + this->Electron_Pressure_Gradient[c+1]/this->Plasma_Density_m3[c+1]/1.602176634e-19);

            int1 += 0.5 * Dz * (int1_1 + int1_2);

            //the area is only correct within the channel, need to add plume area model
            int2_1 = 1/(enemu1*Simulation_Parameters.Channel_Area_m2);
            int2_2 = 1/(enemu2*Simulation_Parameters.Channel_Area_m2);

            int2 += 0.5 * Dz * (int2_1 + int2_2);

        }

        this->Id = (Simulation_Parameters.Discharge_Voltage_V + int1)/int2;
    }

    void Electrons_Object::Compute_Electric_Field(HypiC::Options_Object Simulation_Parameters){
        for(size_t i=0; i<Simulation_Parameters.nCells; ++i){
            double E = ((this->Id / Simulation_Parameters.Channel_Area_m2 - this->Ion_Current_Density[i])/1.602176634e-19/this->Electron_Mobility[i]/this->Plasma_Density_m3[i] - this->Electron_Pressure_Gradient[i]/1.602176634e-19/this->Plasma_Density_m3[i]);
            this->Electric_Field_V_m[i] = E;
        }
    }

    void Electrons_Object::Solve_Potential(HypiC::Options_Object Simulation_Parameters){
        this->Potential[0] = Simulation_Parameters.Discharge_Voltage_V;

        for(size_t i=1; i<Simulation_Parameters.nCells; ++i){
            double dx = this->Cell_Center[i] - this->Cell_Center[i-1];
            this->Potential[i] = this->Potential[i-1] - 0.5 * dx * (this->Electric_Field_V_m[i] + this->Electric_Field_V_m[i-1]);
        }
    }

    void Electrons_Object::Update_Electron_Energy(HypiC::Options_Object Simulation_Parameters, HypiC::Rate_Table_Object Loss_Rates){
        std::vector<double> diag(Simulation_Parameters.nCells,1);
        std::vector<double> diag_low(Simulation_Parameters.nCells-1,1);
        std::vector<double> diag_up(Simulation_Parameters.nCells-1,1);
        std::vector<double> B(Simulation_Parameters.nCells,1);
        std::vector<double> Energy_new;
        double je_sheath;
        double T0;

        //setting boundary conditions
        diag[0] = 1;
        diag_up[0] = 0;
        diag[Simulation_Parameters.nCells-1] = 1;
        diag_low[Simulation_Parameters.nCells-2] = 0;

        //anode neuman for sheath boundary        
        /*B[0] = 0;
        diag[0] = 1/this->EnergyDensity[0];
        diag_up[0] = -1/this->EnergyDensity[1];*/

        //anode boundary condition 
        B[0] = 1.5 * Simulation_Parameters.Initial_Anode_Temperature_eV * this->Plasma_Density_m3[0];
        //cathode boundary condition
        B[Simulation_Parameters.nCells-1] = 1.5 * Simulation_Parameters.Initial_Cathode_Temperature_eV * this->Plasma_Density_m3[Simulation_Parameters.nCells-1];

        //calculate source terms
        //there should be a timestep term in here too 
        std::vector<double> OhmicHeating(Simulation_Parameters.nCells,0);
        std::vector<double> Collisional_Loss(Simulation_Parameters.nCells,0);
        std::vector<double> WallPowerLoss(Simulation_Parameters.nCells,0);
        for(size_t i=1; i<Simulation_Parameters.nCells-1; ++i){ //set limits to disclude anode and cathode  
            OhmicHeating[i] = this->Plasma_Density_m3[i] * this->Electron_Velocity_m_s[i] * this->Electric_Field_V_m[i];
            Collisional_Loss[i] = this->Plasma_Density_m3[i] * this->Neutral_Density_m3[i] * Loss_Rates.interpolate(1.5 * this->Electron_Temperature_eV[i]);
            
            if (this->Cell_Center[i] <= Simulation_Parameters.Channel_Length_m){
                WallPowerLoss[i] = 7.5e6 * this->Electron_Temperature_eV[i] * exp( - 40 / (3 *this->Electron_Temperature_eV[i]));
            } else{
                WallPowerLoss[i] = 1.5e7 * this->Electron_Temperature_eV[i] * exp( - 40 / (3 *this->Electron_Temperature_eV[i]));
            }
            B[i] = Simulation_Parameters.dt * (OhmicHeating[i] - this->Plasma_Density_m3[i] * WallPowerLoss[i] - Collisional_Loss[i]) + this->EnergyDensity[i]; 
        }
        
        //calculate fluxes (setting up linear matrix)
        /*we have two terms we need to care about here are the convective flux (div 5/3 * flux*energy)
        and the heat flux (div kappa nabla energy). By applying the divergence theorem we can remove the div to 
        make the equation only based on the fluxes at the edges

        For the convective flux, use upwinding, which means that we use the value in the same direction as the velocity
        for 1D, if velocity is positive, use the cell to the left, if negative, use the cell to the right
        do this for both edges of each cell. Then just need a factor of 5/3*ue for the respective cell's energy density

        For the heat flux, use upwinding for which cell's kappa to apply and then use a central difference scheme to finite
        difference the gradient in the energy density which adds a geometric factor. 
        */

        //loop over main body
        for(size_t i=1; i<Simulation_Parameters.nCells - 1; ++i){
            if (this->Electron_Velocity_m_s[i] > 0){
                //upwind from 0 to 1 and 1 to 2
                diag_low[i-1] = (-5.0/3.0) * this->Electron_Velocity_m_s[i-1] - this->Electron_Thermal_Conductivity[i-1] / this->Grid_Step;
                diag[i] = 1 + (5.0/3.0) * this->Electron_Velocity_m_s[i] + (this->Electron_Thermal_Conductivity[i-1] + this->Electron_Thermal_Conductivity[i]) / this->Grid_Step;
                diag_up[i] = -this->Electron_Thermal_Conductivity[i] / this->Grid_Step;
            } else{
                //upwind from 1 to 0 and 2 to 1
                diag_low[i-1] =  -1.0 * this->Electron_Thermal_Conductivity[i] / this->Grid_Step;
                diag[i] = 1 + (-5.0/3.0) * this->Electron_Velocity_m_s[i] + (this->Electron_Thermal_Conductivity[i] + this->Electron_Thermal_Conductivity[i+1] ) / this->Grid_Step;
                diag_up[i] = (5.0/3.0) * this->Electron_Velocity_m_s[i+1] -this->Electron_Thermal_Conductivity[i+1] / this->Grid_Step;
            }
            diag_low[i-1] *= Simulation_Parameters.dt / this->Grid_Step;
            diag[i] *= Simulation_Parameters.dt / this->Grid_Step;
            diag_up[i] *= Simulation_Parameters.dt / this->Grid_Step;
            

        }

        //apply sheath boundary condition 
        //neglect heat flux at left edge
        /*if (this->Electron_Velocity_m_s[0] > 0){ 
            //upwind from 1 to 2
            diag[0] = (5.0/3.0) * this->Electron_Velocity_m_s[0] + this->Electron_Thermal_Conductivity[0] / this->Grid_Step;
            diag_up[0] = -this->Electron_Thermal_Conductivity[0] / this->Grid_Step;
        } else{
            //upwind from 2 to 1
            diag[0] = this->Electron_Thermal_Conductivity[1] / this->Grid_Step;
            diag_up[0] = (5.0/3.0) * this->Electron_Velocity_m_s[1] - this->Electron_Thermal_Conductivity[1] / this->Grid_Step;
        }

        //handle the left edge
        T0 = (2.0/3.0) * this->EnergyDensity[0] / this->Plasma_Density_m3[0];
        je_sheath = (this->Id / Simulation_Parameters.Channel_Area_m2) - this->Ion_Current_Density[0];
        //the back term (1-log) of this is the sheath potential 
        diag[0] = (4.0/3.0) * je_sheath / (-1.602176634e-19 * this->Plasma_Density_m3[0]) * (1.0 - log(std::min(1.0, je_sheath / (1.602176634e-19 * this->Plasma_Density_m3[0] * sqrt(8 * 1.602176634e-19 * T0/ (M_PI * 9.10938356e-31)) / 4))));
        */


        //call matrix solver, might want to use Thomas https://www.quantstart.com/articles/Tridiagonal-Matrix-Algorithm-Thomas-Algorithm-in-C/
        Energy_new = HypiC::Thomas_Algorithm(diag_low, diag, diag_up, B);

        //limit to a minimum temperature and assign
        for(size_t i=0; i<Simulation_Parameters.nCells; ++i){
            this->EnergyDensity[i] = std::max(Energy_new[i], 1.5 * this->Plasma_Density_m3[i] * Simulation_Parameters.Min_Electron_Temperature_eV);
        }
    }
}