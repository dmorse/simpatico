#ifndef HYBRID_MD_MOVE_H
#define HYBRID_MD_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>  // base class
#include <util/containers/DArray.h>   // member template
#include <util/space/Vector.h>         // member template parameter

namespace McMd
{

   using namespace Util;

   class McSystem;
   class MdSystem;
   

   /**
   * HybridMdMove is a hybrid Molecular Dynamics MC move.
   *
   * \ingroup McMove_Module MD_Module
   */
   class HybridMdMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      *
      * Constructs a component MdSystem object.
      */
      HybridMdMove(McSystem& system);
   
      /**
      * Destructor.
      */
      ~HybridMdMove();
   
      /**
      * Read nStep, dt, skin, maxNPair from file.
      */
      virtual void readParam(std::istream& in);
   
      /**
      * Generate, attempt and accept or reject a move.
      */
      bool move();
   
   private:
  
      /// MdSystem object used for MD integration
      MdSystem      *mdSystemPtr_;  
   
      /// Array to store old atomic positions.
      DArray<Vector> oldPositions_; 

      /// Number of Md steps per Hybrid MD move
      int            nStep_;
   };

}      
#endif
