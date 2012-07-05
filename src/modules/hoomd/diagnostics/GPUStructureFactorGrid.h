#ifndef GPU_STRUCTURE_FACTOR_GRID_H
#define GPU_STRUCTURE_FACTOR_GRID_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/system/StructureFactorGrid.h>

/** 
 * This class calculates the static structure factor on the GPU.
 */
namespace McMd
{

   using namespace Util;

   class GPUStructureFactorGrid : public StructureFactorGrid
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      GPUStructureFactorGrid(System &system);

      /**
      * Add particles to StructureFactor accumulators.
      *
      * \param iStep step counter
      */
      void sample(long iStep);
    
   };

}
#endif
