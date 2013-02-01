#ifndef MCMD_CFB_HOMO_REPTATION_MOVE_H
#define MCMD_CFB_HOMO_REPTATION_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/CfbEndBase.h>   // base class

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * Configuration bias reptation move for a Linear species.
   *
   * A reptation move uses a configuration bias algorithm to add a segment 
   * to one end of a chain (the head) and delete one from the other (the 
   * tail).  The direction of the move (i.e., which end is the head) is 
   * chosen at random. The only parameter is nTrial, the number of trial 
   * positions for the new position of the head monomer.
   *
   * \ingroup McMd_McMove_Module
   */
   class CfbHomoReptationMove : public CfbEndBase
   {
   
   public:
   
      /**
      * Constructor. 
      */
      CfbHomoReptationMove(McSystem& system);
   
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
  
      #if 0 
      /**
      * Return the associated junction factor, if present. Introduced
      * diblock copolymers.
      */
      double junctionFactor(int headMol, int tailId, int endSign);
      #endif
   
   };

}      
#endif
