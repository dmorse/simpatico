#ifndef TOOLS_MOLECULE_H
#define TOOLS_MOLECULE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>     // members
#include <util/global.h>           // error handling

namespace Tools
{

   using namespace Util;

   struct Atom;
   class Species;

   /**
   * An Molecule has a sequence of atoms, and belongs to an Species.
   *
   * The Molecule class provides read-only access atoms within a 
   * molecule, but does not provide methods to modify its own state. 
   * Molecule relies on Species, a friend class, to control its
   * state. 
   * 
   * \ingroup Tools_Chemistry_Module
   */
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
      * Return Atom identified by index.
      * 
      * \param id atom index within the molecule
      */
      Atom& atom(int id) const;

      /**
      * Is this molecule active?
      */
      bool isActive() const;

   private:

      /**
      * Pointer to an array of pointers to atoms.
      *
      * This member points to a subblock of a larger Atom* array
      * owned by the parent Species.
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

   // Inline functions

   inline Species& Molecule::species() const
   {
      assert(speciesPtr_);
      return *speciesPtr_;
   }

   inline int Molecule::id() const
   {  return id_; }

   inline Atom& Molecule::atom(int id) const
   {
      assert(id >= 0);
      assert(id < nAtom_);
      assert(atoms_[id]);
      return *(atoms_[id]);
   }

   inline bool Molecule::isActive() const
   {  return (bool)nAtom_; }

}
#endif
