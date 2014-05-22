#ifndef MDPP_MOLECULE_H
#define MDPP_MOLECULE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

//#define UTIL_32BIT

#include <util/space/Vector.h>     // members
#include <util/global.h>           // error handling

namespace MdPp
{

   using namespace Util;

   class Atom;
   class Species;

   class Molecule
   {
   public:

      /**
      * Default constructor.
      */
      Molecule();

      /**
      * Return parent Species by reference.
      */ 
      Species& species() const;
   
      /**
      * Return id of this Molecule within its species.
      */ 
      int id() const;
      
      /**
      * Return Atom number id by reference.
      */ 
      Atom& atom(int id) const;
   
   private:
  
      /**
      * Pointer to an array of pointers to atoms.
      */ 
      Atom** atoms_;

      /**
      * Pointer to parent Species.
      */
      Species* speciesPtr_;

      /**
      * Index of molecule within its Species.
      */
      int id_;
   
      /**
      * Number of atoms within this molecule.
      */
      int nAtom_;

   //friends:
   
      friend class Species;
   
   };

   inline Species& Molecule::species() const
   { 
      assert(speciesPtr_);
      assert(id < nAtom_);
      return *speciesPtr_;
   }

   inline int Molecule::id() const
   {  return id_; }
   
   inline Atom& Molecule::atom(int id) const
   { 
      assert(id >=0);
      assert(id < nAtom_);
      return *(atoms_[id]);
   }
   
}
#endif
