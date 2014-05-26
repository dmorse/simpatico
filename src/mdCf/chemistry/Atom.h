#ifndef MDPP_ATOM_H
#define MDPP_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

//#define UTIL_32BIT

#include <util/space/Vector.h>            // members

namespace MdCf
{

   using namespace Util;

   /**
   * A point particle in an MD simulation.
   */
   struct Atom
   {

   public:

      Vector position;
      Vector velocity;
      int typeId;
      int id;
      int speciesId;
      int moleculeId;
      int atomId;

   };

}
#endif
