#pragma once 
#include <memory>
#include <string>
#include <vector>


namespace HypiC{

    //template <class fp_type>
    class Rate_Table_Object{
        private:
            
        public:
            size_t _nEntries = 0;
            std::vector<double> _Energies;
            std::vector<double> _Rates;


            Rate_Table_Object();
            virtual ~Rate_Table_Object();
            void Read_Table(std::string filename);
            double interpolate(double TeV);    
    };
}