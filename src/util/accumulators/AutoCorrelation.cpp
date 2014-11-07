#ifndef UTIL_AUTOCORRELATION_CPP
#define UTIL_AUTOCORRELATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrelation.tpp"  
#include "AutoCorrStage.tpp"  
#include <util/space/Vector.h>
#include <util/space/Tensor.h>

#include <complex>

#include <string>

namespace Util
{
   template class AutoCorrelation<double, double>;
   template class AutoCorrelation<std::complex<double>, std::complex<double> >;
   template class AutoCorrelation<Vector, double>;
   template class AutoCorrelation<Tensor, double>;
   template class AutoCorrStage<double, double>;
}
#endif
