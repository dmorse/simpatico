#ifndef MCMD_RIGID_DISPLACE_MOVE_H
#define MCMD_RIGID_DISPLACE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>  // base class
#include <util/containers/DArray.h>   // member
#include <util/space/Vector.h>         // member template parameter
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * Random rigid displacement of a molecule.
   *
   * \ingroup McMove_Module MD_Module
   */
   class RigidDisplaceMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      */
      RigidDisplaceMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParam(std::istream& in);
   
      /**
      * Generate, attempt and accept or reject a move.
      */
      virtual bool move();

   private:

      /// Array of old positions.
      DArray<Vector> oldPositions_;

      /// Maximum magnitude of displacement
      double delta_;

      /// Integer Id of Species.
      int    speciesId_;
   
      /// Number of atoms per molecule for this species.
      int    nAtom_;
   
   };

}      
#endif
