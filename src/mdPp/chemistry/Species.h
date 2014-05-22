#ifndef MDPP_SPECIES_H
#define MDPP_SPECIES_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Molecule.h"
#include <util/containers/DArray.h>
#include <util/containers/ArrayIterator.h>

namespace MdPp
{

   class Atom;
   using namespace Util;

   class Species {
   public:
   
     typedef ArrayIterator<Molecule> MoleculeIterator;

     Species();
   
     void setId(int id);
   
     void initialize();
   
     void addAtom(Atom& atom);

     void clear();
   
     void begin(MoleculeIterator& iterator);
   
     Molecule& molecule(int i);
   
     /// Return integer id for this species.
     int id() const;
   
     /// Number of molecules.
     int size() const;
   
     void isValid();
   
   private:
   
     DArray<Atom*> atomPtrs_;
     DArray<Molecule> molecules_;
   
     /// Species index
     int id_;
    
     /// Number of atoms per molecule
     int nAtom_;
   
     /// Maximum number of molecules.
     int capacity_;
   
     /// Actual number of molecules = maximum molecule id + 1
     int size_;
   
   //friends:
   
     friend std::istream& operator >> (std::istream& in, Species& species);
     friend std::ostream& operator << (std::ostream& out, Species& species);
   
   };

   /// Return a specific molecule. 
   inline Molecule& Species::molecule(int i)
   {
      assert(i >=0);
      assert(i < size_);
      return molecules_[i]; 
   }
 
   /// Return integer id for this species.
   inline int Species::id() const
   {  return id_; }
 
   /// Number of molecules.
   inline int Species::size() const
   {  return size_; }
 
}
#endif
