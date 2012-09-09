#ifndef MCMD_HOMOPOLYMER_SEMI_GRAND_MOVE_H
#define MCMD_HOMOPOLYMER_SEMI_GRAND_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>   // base class

namespace McMd
{

   using namespace Util;

   class HomopolymerSG;
   class McSystem;

   /**
   * A move that changes the type of a HomopolymerSG molecule.
   *
   * \ingroup McMd_McMove_Module
   */
   class HomopolymerSemiGrandMove : public SystemMove
   {
   
   public:
   
      /**
      * Constructor. 
      */
      HomopolymerSemiGrandMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
   
   protected:
   
      /// Integer index for molecular species.
      int speciesId_;

      /// Pointer to instance of HomopolymerSG.
      HomopolymerSG* speciesPtr_;

   };

}      
#endif
