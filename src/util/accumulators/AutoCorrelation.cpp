#ifndef UTIL_AUTOCORRELATION_CPP
#define UTIL_AUTOCORRELATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrelation.h"  
#include "AutoCorrStage.tpp"  

#include <string>

namespace Util
{

   template class AutoCorrelation<double, double>;
   template class AutoCorrStage<double, double>;
}
#endif
