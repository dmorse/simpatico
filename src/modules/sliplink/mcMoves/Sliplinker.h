#ifndef SLIPLINKER_H
#define SLIPLINKER_H

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
   * Move to create and destroy slip-springs.
   * 
   * \ingroup McMove_Module MD_Module
   */
   class Sliplinker : public SystemMove
   {

   public:
   
      /// Constructor.
      Sliplinker(McSystem& system);

      /// Read cutoff and probability.
      virtual void readParameters(std::istream& in);
 
      /// Create or destroy slip-springs.
      virtual bool move();

   private:

      /// Array to hold neighbors returned by a CellList.
      mutable CellList::NeighborArray neighbors_;
      
      double cutoff_;
      double mu_;
      int speciesId_;
   
   };

}

#endif
