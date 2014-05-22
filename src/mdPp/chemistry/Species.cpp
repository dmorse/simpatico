#ifndef MDPP_SPECIES_CPP
#define MDPP_SPECIES_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"

namespace MdPp
{

   using namespace Util;
   
   void Species::addAtom(Atom& atom) 
   {
   
      // Check that atom.speciesId = id_;
      // Check that moleculeId >= 0 and < capacity_
      // Check that atomId >=0 and < nAtom_
      // int i = (moleculeId*nAtom_) + atomId_;
      // atomPtrs_[i] = &atom;
      // if (moleculeId > size_) {
      //    size_ = moleculeId;
      //}
   }
   
   void Species::isValid() 
   {
      // Check that atom pointers have been set for all atoms in all
      // molecules with moleculeId < size_, and that AtomContext info
      // in each atom is consistent with data in Species and Molecule.
   }

}
#endif
