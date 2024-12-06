#include <cmath>
#include <fstream>
#include "HypiCpp.hpp"
#include <iostream>
#include <random>

namespace HypiC{
    //class functions
    Options_Object::Options_Object()
    {
        
    }
    
    // Destructor
    //template <class fp_type>
    //Particles_Object<fp_type>::~Particles_Object()
    Options_Object::~Options_Object()
    {
        
    }
    void Options_Object::Read_Input(std::string Filename){
        std::string line;
        size_t pos;
    
        //open the file and read line by line
        //open the file 
        std::ifstream f(Filename);
        
        //grab line, allow for comments in the input file, each line should end with //
        std::getline(f, line);
        pos = line.find_first_of("//");
        //line 1 is nIterations
        this->nIterations = std::stoul(line.substr(0, pos));
        //rinse and repeat, next line is Output_Interval
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Output_Interval = std::stoul(line.substr(0, pos));
        //dt
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->dt = std::stod(line.substr(0, pos));
        //nCells
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->nCells = std::stoul(line.substr(0, pos));
        //Domain Length
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Domain_Length_m = std::stod(line.substr(0, pos));
        //Channel Length
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Channel_Length_m = std::stod(line.substr(0, pos));
        //Channel Area
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Channel_Area_m2 = std::stod(line.substr(0, pos));
        //Discharge Voltage
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Discharge_Voltage_V = std::stod(line.substr(0, pos));
        //Mass Flow Rate
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Mass_Flow_Rate_kg_s = std::stod(line.substr(0, pos));
        //Number of neutral particles
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->N_Neutrals = std::stoul(line.substr(0, pos));
        //Number of ion particles
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->N_Ions = std::stoul(line.substr(0, pos));
        //Neutral Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Neutral_Temperature_K = std::stod(line.substr(0, pos));
        //Ion Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Ion_Temperature_K = std::stod(line.substr(0, pos));
        //Electron Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Max_Electron_Temperature_eV = std::stod(line.substr(0, pos));
        //Anode Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Anode_Temperature_eV = std::stod(line.substr(0, pos));
        //Cathode Temperature
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Cathode_Temperature_eV = std::stod(line.substr(0, pos));
        //Initial_Min_Ion_Denisty
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Min_Ion_Density = std::stod(line.substr(0, pos));
        //Max Ion Density
        std::getline(f, line);
        pos = line.find_first_of("//");
        this->Initial_Max_Ion_Density = std::stod(line.substr(0, pos));
    };
    
    double Maxwellian_Sampler(double mu, double sigma){
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution Dist{mu, sigma};
        //sample from a Maxwellian 
        return Dist(gen);
    };

    double Initial_Magnetic_Field(double B_max, double Lch, double z){
        double B;
        if (z < Lch){
            B = B_max * exp(-0.5 * pow(((z-Lch)/0.011),2));
        } else {
            B = B_max * exp(-0.5 * pow(((z-Lch)/0.018),2));
        }
        return B;
    }

    double Initial_Electron_Density(double z, double n_min, double n_max, double Vd, double mdot, double Lch){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/
        return  sqrt(Vd/300) * (mdot / 5e-6) * (n_min + (n_max - n_min)*exp(-1.0*pow(((3.0*z/Lch) - 1.5),2)));
    };

    double Initial_Neutral_Density(double z, double mass, double un, double mdot, double Lch, double Ach){
        double n_anode;
        double n_cathode;

        n_anode = mdot / (mass * un * Ach);
        n_cathode = 0.01 * n_anode;
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/
        return 0.5 * (n_anode + n_cathode + (n_anode - n_cathode) * tanh((z - 0.5 * Lch)/(Lch / 6)));
    };

    double Initial_Electron_Temperature(double z, double Te_Anode, double Te_Cathode, double Te_Max, double Lch, double z_max){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/
        double Te_min = std::min(Te_Anode, Te_Cathode);
        return (1 - (z/z_max)) * Te_Anode + (z/z_max) * Te_Cathode + (Te_Max - Te_min) * exp(-1.0*pow(((3.0*z/Lch) - 3),2));
    };

    double Initial_Ion_Bulk_Velocity(double Te_Anode, double Vd, double z, double Lch, double z_max){
        //from HallThruster.jl 
        //see https://um-pepl.github.io/HallThruster.jl/dev/initialization/

        double u_Bohm = -1 * sqrt(1.6e-19 * Te_Anode / (1.6e-27 * 131.29));//e=1.6e-19 C, mi_Xe = 1.6e-27 (mp) * 131.29 (amu Xe)
        double u_max = sqrt(2 * 1.6e-19 * Vd / (1.6e-27 * 131.29));//e=1.6e-19 C, mi_Xe = 1.6e-27 (mp) * 131.29 (amu Xe)
        
        if (z<Lch){
            return u_Bohm + (2/3*(u_max - u_Bohm))*pow(z/Lch,2);
        } else{
            return (2* u_max / 3 + u_Bohm / 3)*(1-(z-Lch)/(z_max-Lch)) + u_max * ((z-Lch)/(z_max-Lch));
        } 
    };

    //for the initialization of particles, I think we can physically distribute the particles evenly
    //velocities come from maxwellian
    //weights come from density 
    //electron density and electron temperature come from the inital distributions
    //electric field can be set to 0 then updated (update electrons first?).

    HypiC::Particles_Object Initialize_Neutrals(HypiC::Options_Object Inputs)
    {
        int Particles_per_cell;
        double z;
        double z_particle;
        double dz;
        double nn;
        double ne;
        double Te;
        double v;
        double w;
        double mass;
        double kb;
        double un;

        //calculate grid increment (assuming uniformly spaced)
        dz = Inputs.Domain_Length_m / Inputs.nCells;

        //create the particle class instance
        HypiC::Particles_Object Neutrals = HypiC::Particles_Object();
        Neutrals._IonizationDirection = -1.0;//neutrals are removed due to ionization
        mass = 131.29 * 1.66053907e-27;//for Xe
        kb = 1.380649e-23;
        srand(time(NULL));

        //loop over grid cells
        for(size_t c=0; c<Inputs.nCells; ++c){
            //calculate number of particles per cell and the position
            Particles_per_cell = Inputs.N_Neutrals / Inputs.nCells;
            z = dz * c;
            //loop over the particles in the cell
            for(size_t p=0; p<Particles_per_cell; ++p){
                //uniformly sample across cell. 
                z_particle = z + dz * (rand()/RAND_MAX);
                ne = HypiC::Initial_Electron_Density(z_particle, Inputs.Initial_Min_Ion_Density,
                Inputs.Initial_Max_Ion_Density, Inputs.Discharge_Voltage_V, Inputs.Mass_Flow_Rate_kg_s, 
                Inputs.Channel_Length_m);
                Te = HypiC::Initial_Electron_Temperature(z_particle, Inputs.Initial_Anode_Temperature_eV,
                Inputs.Initial_Cathode_Temperature_eV, Inputs.Initial_Max_Electron_Temperature_eV, 
                Inputs.Channel_Length_m, Inputs.Domain_Length_m);
                //calculate initial velocity
                //sample from maxwellian
                un = sqrt(2 * kb * Inputs.Initial_Neutral_Temperature_K / (M_PI * mass));
                v = HypiC::Maxwellian_Sampler(un,sqrt(kb * Inputs.Initial_Neutral_Temperature_K / mass));
                //calculate initial weight
                //see https://smileipic.github.io/Smilei/Understand/algorithms.html
                //use neutral density profile
                nn = Initial_Neutral_Density(z_particle, mass, un, Inputs.Mass_Flow_Rate_kg_s, Inputs.Channel_Length_m,
                Inputs.Channel_Area_m2);
                w = (nn / Particles_per_cell) * dz;

                //Add the particle
                Neutrals.Add_Particle(z_particle, v, w, 0, ne, Te);
            }
        }

        //return
        return Neutrals;
    };

    HypiC::Particles_Object Initialize_Ions(HypiC::Options_Object Inputs)
    {
        int Particles_per_cell;
        double z;
        double z_particle;
        double dz;
        double ne;
        double Te;
        double v;
        double w;
        double mass;
        double kb;

        //calculate grid increment (assuming uniformly spaced)
        dz = Inputs.Domain_Length_m / Inputs.nCells;

        //create the particle class instance
        HypiC::Particles_Object Ions = HypiC::Particles_Object();
        Ions._IonizationDirection = 1.0;//ions are added due to ionization
        mass = 131.29 * 1.66053907e-27;//for Xe
        kb = 1.380649e-23;
        Ions._ChargetoMassRatio = 1.602176634e-19 / mass; 
        srand(time(NULL));

        //loop over grid cells
        for(size_t c=0; c<Inputs.nCells; ++c){
            //calculate number of particles per cell and the position
            Particles_per_cell = Inputs.N_Ions / Inputs.nCells;
            z = dz * c;
            
            //loop over the particles in the cell
            for(size_t p=0; p<Particles_per_cell; ++p){
                //uniformly sample
                //later add to the position 
                z_particle = z + dz * (rand()/RAND_MAX);

                ne = HypiC::Initial_Electron_Density(z_particle, Inputs.Initial_Min_Ion_Density,
                Inputs.Initial_Max_Ion_Density, Inputs.Discharge_Voltage_V, Inputs.Mass_Flow_Rate_kg_s, 
                Inputs.Channel_Length_m);
                Te = HypiC::Initial_Electron_Temperature(z_particle, Inputs.Initial_Anode_Temperature_eV,
                Inputs.Initial_Cathode_Temperature_eV, Inputs.Initial_Max_Electron_Temperature_eV, 
                Inputs.Channel_Length_m, Inputs.Domain_Length_m);
                //calculate initial velocity
                //sample from a maxwellian
                v = HypiC::Maxwellian_Sampler(HypiC::Initial_Ion_Bulk_Velocity(Inputs.Initial_Anode_Temperature_eV, Inputs.Discharge_Voltage_V,
                z_particle, Inputs.Channel_Length_m, Inputs.Domain_Length_m), sqrt(kb * Inputs.Initial_Ion_Temperature_K / mass));
                if (v>3e8){
                    std::cout << z_particle << "\n";
                    std::cout << HypiC::Initial_Ion_Bulk_Velocity(Inputs.Initial_Anode_Temperature_eV, Inputs.Discharge_Voltage_V,
                    z_particle, Inputs.Channel_Length_m, Inputs.Domain_Length_m) << "\n";
                    
                }
                //calculate initial weight, for singly charged only can use the electron number density
                //see https://smileipic.github.io/Smilei/Understand/algorithms.html
                w = (ne / Particles_per_cell) * dz;

                //Add the particle
                Ions.Add_Particle(z, v, w, 0, ne, Te);
            }
        }

        //return
        return Ions;
    };

    HypiC::Electrons_Object Initialize_Electrons(HypiC::Options_Object Inputs)
    {
        int Particles_per_cell;
        double z;
        double z_particle;
        double dz;
        double ne;
        double Te;
        double v;
        double w;
        double mass;
        double kb;
        double EnergyDensity;
        double B;
        double f;
        double Efield;

        //calculate grid increment (assuming uniformly spaced)
        dz = Inputs.Domain_Length_m / Inputs.nCells;

        //create the particle class instance
        HypiC::Electrons_Object Electrons = HypiC::Electrons_Object();
        Electrons.Grid_Step = dz;
        mass = 131.29 * 1.66053907e-27;//for Xe
        kb = 1.380649e-23;
        srand(time(NULL));

        //loop over grid cells
        for(size_t c=0; c<Inputs.nCells; ++c){
            //calculate number of particles per cell and the position
            z = 0.5 * dz + dz * c;

            ne = HypiC::Initial_Electron_Density(z, Inputs.Initial_Min_Ion_Density,
            Inputs.Initial_Max_Ion_Density, Inputs.Discharge_Voltage_V, Inputs.Mass_Flow_Rate_kg_s, 
            Inputs.Channel_Length_m);
            Te = HypiC::Initial_Electron_Temperature(z, Inputs.Initial_Anode_Temperature_eV,
            Inputs.Initial_Cathode_Temperature_eV, Inputs.Initial_Max_Electron_Temperature_eV, 
            Inputs.Channel_Length_m, Inputs.Domain_Length_m);
            B = HypiC::Initial_Magnetic_Field(0.015, Inputs.Channel_Length_m, z);

            EnergyDensity = ne * Te * 3/2;

            v = HypiC::Initial_Ion_Bulk_Velocity(Inputs.Initial_Anode_Temperature_eV, Inputs.Discharge_Voltage_V,
            z, Inputs.Channel_Length_m, Inputs.Domain_Length_m);

            // Utilizing Case 2 from 
            // https://0534de96-08f4-4c2b-82e3-b2a3f3216551.filesusr.com/ugd/8243e7_030af6befce7412ca1232089b304903e.pdf
            if (z <= Inputs.Channel_Length_m) {
                f = 0.5e7;
            } else {
                f = 1e7;
            }

            Efield = 2;

            //Add the particle
            Electrons.Add_Electron(ne, Te, B, EnergyDensity, v, f, Efield, z);
        }

        //return
        return Electrons;
    };

}