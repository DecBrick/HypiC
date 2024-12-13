#include "Rate_Table.hpp"
#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>


namespace HypiC{

    // Default constructor
    //template <class fp_type>
    //Rate_Table_Object<fp_type>::Rate_Table_Object()
    Rate_Table_Object::Rate_Table_Object()
    {
        //std::cout << "Hello from Rate_Table_Constructor!\n";
    }
    
    // Destructor
    //template <class fp_type>
    //Rate_Table_Object<fp_type>::~Rate_Table_Object()
    Rate_Table_Object::~Rate_Table_Object()
    {
        //delete[] this->_Energies;
        //delete[] this->_Rates;
        //std::cout << "Byebye from Rate_Table_Destructor~\n";
    }


    //Read Table file
    //template <class fp_type>
    //void Rate_Table_Object<fp_type>::Read_Table(std::string filename){
    void Rate_Table_Object::Read_Table(std::string filename){
        size_t count = 0;
        //open the file 
        std::ifstream f(filename);
        //loop over each line
        for (std::string line; std::getline(f, line);){
            //find the first space/tab
            size_t pos = line.find_first_of(" \t");
            //interpret the first part of the line as the energy 
            this->_Energies.push_back(std::stod(line.substr(0, pos)));
            //interpret the last part of the line as the rate 
            this->_Rates.push_back(std::stod(line.substr(pos + 1)));
            //count 
            count+=1;
        }
        //save the number of entires
        this->_nEntries = count; 

    }


    //Interpolation method
    
    //template <class fp_type>
    //fp_type Rate_Table_Object<fp_type>::interpolate(fp_type TeV){
    double Rate_Table_Object::interpolate(double TeV){
        
        double m;

        //check bounds and extrapolate.
        if (TeV <= this->_Energies[0]){
            return this->_Rates[0];
        }
        if (TeV >= this->_Energies[this->_nEntries-1]){//-1 is to account for the 0 based-indexing
            return this->_Rates[this->_nEntries-1];
        }

        //compare index by index to find left (assumes sorted)
        for (size_t i = 1; i < this->_nEntries; ++i){
            if (TeV <= this->_Energies[i]){
                m = (this->_Rates[i] - this->_Rates[i-1]) / (this->_Energies[i] - this->_Energies[i-1]);
                return m * (TeV - this->_Energies[i-1]) + this->_Rates[i-1];
            }
        }
        //default return statement to make compiler happy
        return -1;
    }
    
}