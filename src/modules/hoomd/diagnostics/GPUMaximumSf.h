#ifndef GPU_MAXIMUM_SF_H
#define GPU_MAXIMUM_SF_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/system/MaximumSf.h>
#include "GPUStructureFactorGrid.cuh"

/** 
 * This class calculates the static structure factor on the GPU 
 * and finds the maximum structure factor on CPU.
 */
namespace McMd
{

   using namespace Util;

   class GPUMaximumSf : public MaximumSf
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      GPUMaximumSf(System &system);

      /**
      * Add particles to StructureFactor accumulators
      * and evaluate maximum StructureFactor at given iStep
      *
      * \param iStep step counter
      */
      void sample(long iStep);
    
   };

}
#endif
