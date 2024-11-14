#pragma once 
#include <memory>
#include <string>


namespace HypiC{

    //template <class fp_type>
    class Rate_Table_Object{
        private:
            size_t _nEntries = 0;
        public:
            std::unique_ptr< double [] > _Energies = nullptr;
            std::unique_ptr< double [] > _Rates = nullptr;


            Rate_Table_Object();
            virtual ~Rate_Table_Object();
            void Read_Table(std::string filename);
            //fp_type interpolate(fp_type TeV);    
    };
}