#include "Particles.hpp"

namespace HypiC{

    // Default constructor
    template <class fp_type>
    Particles_Object<fp_type>::Particles_Object() :
         HypiC<fp_type>::Particles_Object()
    {
        std::cout << "Hello from Constructor!\n";
    }
    
    // Destructor
    template <class fp_type>
    Particles_Object<fp_type>::~Particles_Object()
    {
        std::cout << "Byebye from  Destructor~\n";
    }


    //Add particle method
    template <class fp_type>
    void Particles_Object<fp_type>::Add_Particle(){
        std::cout << "Hello From Add Particle Method\n";
    }


    //Remove particle method
    template <class fp_type>
    void Particles_Object<fp_type>::Remove_Particle(){
        std::cout << "Hello From Remove Particle Method\n";
    }

    //Update particle method
    template <class fp_type>
    void Particles_Object<fp_type>::Update_Particle(){
        std::cout << "Hello From Update Particle Method\n";
    }
}