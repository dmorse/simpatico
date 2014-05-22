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
  
   // Constructor.
   Species::Species()
    : atomPtrs_(),
      molecules_(),
      id_(-1),
      nAtom_(0),
      capacity_(0),
      size_(0)
   {}
 
   void Species::setId(int id)
   {  id_ = id; }
 
   void Species::initialize()
   {}
 
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
   
   void Species::clear()
   {}
 
   void Species::begin(MoleculeIterator& iterator)
   {  molecules_.begin(iterator); }
 
   void Species::isValid() 
   {
      // Check that atom pointers have been set for all atoms in all
      // molecules with moleculeId < size_, and that AtomContext info
      // in each atom is consistent with data in Species and Molecule.
   }
 
   std::istream& operator >> (std::istream& in, Species& species);
   std::ostream& operator << (std::ostream& out, Species& species);
  
}
#endif
