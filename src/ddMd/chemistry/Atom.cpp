#ifndef DDMD_ATOM_CPP
#define DDMD_ATOM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Atom.h"
#include <ddMd/chemistry/Mask.h>
#include <ddMd/communicate/Plan.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor (private, used by AtomArray).
   */
   Atom::Atom() :
     position_(0.0),
     typeId_(-1),
     localId_(0),
     force_(0.0),
     arrayPtr_(0)
     #ifdef UTIL_32BIT
     , pad_(0)
     #endif
   {}

   /*
   * Assignment (public).
   */
   Atom& Atom::operator= (const Atom& other)
   {
      position_ = other.position_;
      typeId_ = other.typeId_;
      setIsGhost(other.isGhost());
      force_ = other.force_;
      velocity() = other.velocity();
      mask() = other.mask();
      plan() = other.plan();
      setId(other.id());
      return *this;
   }

   /*
   * Reset integer members to null values.
   */
   void Atom::clear()
   {
      typeId_ = -1;
      setIsGhost(false);
      mask().clear();
      plan().clearFlags();
      setId(-1);
   }

}
#endif
