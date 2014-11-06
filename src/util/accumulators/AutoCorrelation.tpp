#ifndef UTIL_AUTOCORRELATION_TPP
#define UTIL_AUTOCORRELATION_TPP

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

   /*
   * Constructor
   *
   * \param blockFactor ratio of block sizes of subsequent stages
   */
   template <typename Data, typename Product>
   AutoCorrelation<Data, Product>::AutoCorrelation()
    : AutoCorrStage()
   {}

   /*
   * Register the creation of a descendant stage.
   */
   template <typename Data, typename Product>
   virtual void 
   AutoCorrelation<Data, Product>::registerDescendant(AutoCorrelation* ptr)
   {}

   /*
   * Output the autocorrelation function
   */
   template <typename Data, typename Product>
   virtual void AutoCorrelation<Data, Product>::output()
   {}

}
#endif
