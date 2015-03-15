#ifndef GC_SLIPLINKMOVE_H
#define GC_SLIPLINKMOVE_H

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>        // base class
#include <mcMd/simulation/System.h>
#include <mcMd/neighbor/CellList.h>
#include <util/global.h>

namespace McMd
{
    
   using namespace Util;
  
   /**
   * Move to evolve the slip-springs.
   * 
   * \ingroup McMove_Module MD_Module
   */
   class GcSliplinkMove : public SystemMove
   {

   public:
   
      /// Constructor.
      GcSliplinkMove(McSystem& system);

      /// Read cutoff and speciesId.
      virtual void readParameters(std::istream& in);
 
      /// Move slip-springs.
      virtual bool move();

   private:

      /// Array to hold neighbors returned by a CellList.
      mutable CellList::NeighborArray neighbors_;
      
      double cutoff_;

      double cutoffSq_;

      double fCreate_;

      double fNotCreate_;

      double mu_;

      int nTrial_;

      int speciesId_;
   
   };

}

#endif
