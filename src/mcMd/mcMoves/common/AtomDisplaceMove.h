#ifndef ATOM_DISPLACE_MOVE_H
#define ATOM_DISPLACE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>        // base class
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * Random displacement of one atom.
   *
   * \ingroup McMove_Module MD_Module
   */
   class AtomDisplaceMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      */
      AtomDisplaceMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParam(std::istream& in);
   
      /**
      * Generate, attempt and accept or reject a move.
      */
      virtual bool move();

   private:

      /// Maximum magnitude of displacement.
      double delta_;

      /// Integer Id of Species.
      int    speciesId_;
   
   };

}      
#endif
