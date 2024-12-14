//define program
#ifndef _HypiC_
#define _HypiC_

//check for OpenMp
#ifdef _OPENMP
    #include <omp.h>
#endif

//add headers for classes 
#include "Particles_Object.hpp"
#include "Electrons_Object.hpp"
#include "Rate_Table.hpp"
#include "Initialize.hpp"
#include "Iteration.hpp"
#include "Input.hpp"
#include "Output.hpp"
#include "math.hpp"
#include "Interpolate.hpp"

#endif // _HypiC_