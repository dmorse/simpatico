#ifndef DDMD_SP_MOLECULE_H
#define DDMD_SP_MOLECULE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>     // members
#include <util/global.h>           // error handling

namespace DdMd
{

   using namespace Util;

   class SpAtom;
   class SpSpecies;

   /**
   * An SpMolecule has a sequence of atoms, and belongs to an SpSpecies.
   *
   * \ingroup DdMd_Sp_Chemistry_Module
   */
   class SpMolecule
   {
   public:

      /**
      * Default constructor.
      */
      SpMolecule();

      /**
      * Return parent SpSpecies by reference.
      */
      SpSpecies& species() const;

      /**
      * Return id of this SpMolecule within its species.
      */
      int id() const;

      /**
      * Return Atom number id by reference.
      */
      SpAtom& atom(int id) const;

   private:

      /**
      * Pointer to an array of pointers to atoms.
      */
      SpAtom** atoms_;

      /**
      * Pointer to parent SpSpecies.
      */
      SpSpecies* speciesPtr_;

      /**
      * Index of molecule within its SpSpecies.
      */
      int id_;

      /**
      * Number of atoms within this molecule.
      */
      int nAtom_;

   //friends:

      friend class SpSpecies;

   };

   // Inline functions

   inline SpSpecies& SpMolecule::species() const
   {
      assert(speciesPtr_);
      assert(id < nAtom_);
      return *speciesPtr_;
   }

   inline int SpMolecule::id() const
   {  return id_; }

   inline SpAtom& SpMolecule::atom(int id) const
   {
      assert(id >=0);
      assert(id < nAtom_);
      return *(atoms_[id]);
   }

}
#endif
