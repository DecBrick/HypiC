#include "Rate_Table.hpp"
#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>


namespace HypiC{

    // Default constructor
    template <class fp_type>
    Rate_Table_Object<fp_type>::Rate_Table_Object() :
        fp_type::Rate_Table_Object()
    {
        std::cout << "Hello from Rate_Table_Constructor!\n";
    }
    
    // Destructor
    template <class fp_type>
    Rate_Table_Object<fp_type>::~Rate_Table_Object()
    {
        delete[] this->_Energies;
        delete[] this->_Rates;
        std::cout << "Byebye from Rate_Table_Destructor~\n";
    }


    //Read Table file
    template <class fp_type>
    void Rate_Table_Object<fp_type>::Read_Table(std::string filename){
        std::cout << "Hello From Read Table Method\n";
        size_t count;
        //open the file 
        std::ifstream f(filename);
        //loop over each line
        for (std::string line; std::getline(f, line);){
            //find the first space
            std::string::size_type pos = line.find(" ");
            //interpret the first part of the line as the energy 
            this->_Energies[count] = std::stod(line.substr(0, pos));
            //interpret the last part of the line as the rate 
            this->_Rates[count] = std::stod(line.substr(pos + 1));
            //count 
            count+=1;
        }
        //save the number of entires
        this->_nEntries = count -1;//-1 to account for final count of loop 
    }


    //Interpolation method
    template <class fp_type>
    fp_type Rate_Table_Object<fp_type>::interpolate(fp_type TeV){
        std::cout << "Hello From interpolate Method\n";
        double m;

        //check bounds and extrapolate.
        if (TeV <= this->_Energies[0]){
            return this->_Rates[0];
        }
        if (TeV >= this->_Energies[this->_nEntries]){
            return this->_Rates[this->_nEntries];
        }

        //compare index by index to find left (assumes sorted)
        for (size_t i = 0; i < this->_nEntries; ++i){
            if (TeV >= this->_Energies[i]){
                m = (this->_Rates[i+1] - this->_Rates[i]) / (this->_Energies[i+1] - this->_Energies[i]);
                return m * (TeV - this->_Energies[i]) + this->_Rates[i];
            }
        }
    }
}