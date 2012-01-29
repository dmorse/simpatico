#ifndef CFB_REBRIDGE_MOVE_H
#define CFB_REBRIDGE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, Jian Qin and David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/CfbRebridgeBase.h>  // base class
#include <util/containers/DArray.h>             // member template
#include <util/space/Vector.h>                   // member template parameter

namespace McMd
{

   using namespace Util;

   class McSystem;
   class Linear;

   /**
   * Config-bias move for internal segment of a flexible linear polymer.
   *
   * \ingroup McMove_Module
   */
   class CfbRebridgeMove : public CfbRebridgeBase
   {
   
   public:
   
      /**
      * Constructor. 
      */
      CfbRebridgeMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParam(std::istream& in);

      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
   
   protected:
   
      /// Array of old positions of temporarily deleted atoms
      DArray<Vector> oldPos_;
    
      /// Integer index for molecular species.
      int speciesId_;
   
      /// Number of particles at end to attempt to regrow
      int nRegrow_;
   
   };

}      
#endif
