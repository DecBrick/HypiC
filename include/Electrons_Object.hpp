#pragma once
#include <memory>//for ptrs 
#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>
#include "Rate_Table.hpp"

namespace HypiC{
    class Electrons_Object{
        private:
            std::vector<double> _ElecDensity;
            std::vector<double> _ChargeDensity;
            size_t _nElectrons = 0;

        public:
            Electrons_Object();
            ~Electrons_Object();

        void Add_Electron(double electron_density, double charge_density);
    };
}
