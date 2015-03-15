#ifndef SLIPLINKMOVE_H
#define SLIPLINKMOVE_H

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
   class SliplinkMove : public SystemMove
   {

   public:

      /// Constructor.
      SliplinkMove(McSystem& system);

      /// Read cutoff and speciesId.
      virtual void readParameters(std::istream& in);

      /// Move slip-springs.
      virtual bool move();

   private:

      /// Array to hold neighbors returned by a CellList.
      mutable CellList::NeighborArray neighbors_;

      double cutoff_;
      int speciesId_;

   };

}

#endif
