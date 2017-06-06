#ifndef MCMD_MOLECULE_H
#define MCMD_MOLECULE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Atom.h"
#ifdef SIMP_BOND
#include "Bond.h"
#endif
#ifdef SIMP_ANGLE
#include "Angle.h"
#endif
#ifdef SIMP_DIHEDRAL
#include "Dihedral.h"
#endif
#include <util/containers/ArrayIterator.h>
#include <util/containers/ConstArrayIterator.h>
#include <util/global.h>

namespace Simp {
   class   Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   class   System;

   /**
   * A physical molecule (a set of covalently bonded Atoms).
   *
   * \ingroup McMd_Chemistry_Module
   */
   class Molecule
   {

   public:

      // Typedefs

      /// Iterator for Atoms within a Molecule.
      typedef ArrayIterator<Atom>  AtomIterator;

      /// Iterator for const Atoms within a Molecule.
      typedef ConstArrayIterator<Atom>  ConstAtomIterator;

      #ifdef SIMP_BOND
      /// Iterator for Bonds within a Molecule.
      typedef ArrayIterator<Bond>  BondIterator;

      /// Iterator for const Bonds within a Molecule.
      typedef ConstArrayIterator<Bond>  ConstBondIterator;
      #endif

      #ifdef SIMP_ANGLE
      /// Iterator for Angles within a Molecule.
      typedef ArrayIterator<Angle>  AngleIterator;

      /// Iterator for const Angles within a Molecule.
      typedef ConstArrayIterator<Angle>  ConstAngleIterator;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Iterator for Dihedrals within a Molecule.
      typedef ArrayIterator<Dihedral>  DihedralIterator;

      /// Iterator for const Dihedrals within a Molecule.
      typedef ConstArrayIterator<Dihedral>  ConstDihedralIterator;
      #endif

      /// Null value for molecule id.
      static const int NullIndex = -1;

      // Member functions

      /**
      * Constructor.
      */
      Molecule();

      /**
      * \name Mutators
      */
      //@{

      /**
      * Set the parent Species.
      *
      * \param species parent Species object.
      */
      void setSpecies(Species &species);

      /**
      * Set the parent System.
      *
      * \param system parent System object.
      */
      void setSystem(System &system);

      /**
      * Set the parent System pointer to null (no System).
      */
      void unsetSystem();

      /**
      * Set the integer index for this molecule.
      *
      * \param id integer array index.
      */
      void setId(int id);

      /**
      * Set the first Atom.
      *
      * \param atom first Atom in molecule.
      */
      void setFirstAtom(Atom &atom);

      /**
      * Set the number of Atoms per molecule.
      *
      * \param nAtom number of Atoms in the molecule.
      */
      void setNAtom(int nAtom);

      #ifdef SIMP_BOND
      /**
      * Set the first Bond.
      *
      * \param bond first Bond in molecule.
      */
      void setFirstBond(Bond &bond);

      /**
      * Set the number of Bonds per molecule.
      *
      * \param nBond number of Bonds in the molecule.
      */
      void setNBond(int nBond);
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Set the first Angle.
      *
      * \param angle first Angle in molecule.
      */
      void setFirstAngle(Angle &angle);

      /**
      * Set the number of Angles per molecule.
      *
      * \param nAngle number of Angles in molecule.
      */
      void setNAngle(int nAngle);
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Set the first Dihedral.
      *
      * \param dihedral first Dihedral in molecule.
      */
      void setFirstDihedral(Dihedral &dihedral);

      /**
      * Set the number of Dihedrals per molecule.
      *
      * \param nDihedral number of Dihedrals in molecule.
      */
      void setNDihedral(int nDihedral);
      #endif

      //@}
      /// \name Accessors
      //@{

      /// Get the molecular Species by reference.
      Species& species() const;

      /// Get the parent System.
      System& system() const;

      /// Get the index for this Molecule (unique in species).
      int id() const;

      /// Get the number of Atoms in this Molecule.
      int nAtom() const;

      #ifdef SIMP_BOND
      /// Get the number of Bonds in this Molecule.
      int nBond() const;
      #endif

      #ifdef SIMP_ANGLE
      /// Get the number of Angles in this Molecule.
      int nAngle() const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Get the number of Dihedrals in this Molecule.
      int nDihedral() const;
      #endif

      /**
      * Get a specific Atom in this Molecule.
      *
      * Returns the atom with local integer index 0 <= localId < nAtom.
      *
      * \param localId local index of desired atom within this Molecule.
      */
      const Atom& atom(int localId) const;

      /**
      * Get a specific Atom in this Molecule.
      *
      * Returns the atom with local integer index 0 <= localId < nAtom.
      *
      * \param localId local index of desired atom within this Molecule.
      */
      Atom& atom(int localId);

      #ifdef SIMP_BOND
      /**
      * Get a specific Bond in this Molecule by non-const reference.
      *
      * Returns the bond with local integer index 0 <= localId < nBond.
      *
      * \param localId local index of desired bond within this Molecule.
      */
      Bond& bond(int localId);

      /**
      * Get a specific Bond in this Molecule by const reference.
      *
      * Returns the bond with local integer index 0 <= localId < nBond.
      *
      * \param localId local index of desired bond within this Molecule.
      */
      const Bond& bond(int localId) const;
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get a specific Angle in this Molecule by non-const reference.
      *
      * Returns the angle with local integer index 0 <= localId < nAngle.
      *
      * \param localId local index of desired angle within this Molecule.
      */
      Angle& angle(int localId);

      /**
      * Get a specific Angle in this Molecule by const reference.
      *
      * Returns the angle with local integer index 0 <= localId < nAngle.
      *
      * \param localId local index of desired angle within this Molecule.
      */
      const Angle& angle(int localId) const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Get a specific Dihedral in this Molecule by reference.
      *
      * Returns the dihedral with local integer index 0 <= localId < nDihedral.
      *
      * \param localId local index of desired dihedral within this Molecule.
      */
      Dihedral& dihedral(int localId);

      /**
      * Get a specific Dihedral in this Molecule by const reference.
      *
      * Returns the dihedral with local integer index 0 <= localId < nDihedral.
      *
      * \param localId local index of desired dihedral within this Molecule.
      */
      const Dihedral& dihedral(int localId) const;
      #endif

      //@}
      /// \name Iterator Interface
      //@{

      /**
      * Set an Molecule::AtomIterator to first Atom in this Molecule.
      *
      * \param iterator on output, points to first Atom.
      */
      void begin(AtomIterator &iterator);

      /**
      * Set an Molecule::ConstAtomIterator to first Atom in this Molecule.
      *
      * \param iterator on output, points to first Atom.
      */
      void begin(ConstAtomIterator &iterator) const;

      #ifdef SIMP_BOND
      /**
      * Set a Molecule::BondIterator to first Bond in this Molecule.
      *
      * \param iterator on output, points to first Bond.
      */
      void begin(BondIterator &iterator);

      /**
      * Set a Molecule::ConstBondIterator to first Bond in this Molecule.
      *
      * \param iterator on output, points to first Bond.
      */
      void begin(ConstBondIterator &iterator) const;
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Set a Molecule::AngleIterator to first Angle in this Molecule.
      *
      * \param iterator on output, points to first Angle.
      */
      void begin(AngleIterator &iterator);

      /**
      * Set a Molecule::ConstAngleIterator to first Angle in this Molecule.
      *
      * \param iterator on output, points to first Angle.
      */
      void begin(ConstAngleIterator &iterator) const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Set a Molecule::DihedralIterator to first Dihedral in this Molecule.
      *
      * \param iterator on output, points to first Dihedral.
      */
      void begin(DihedralIterator &iterator);

      /**
      * Set a Molecule::ConstDihedralIterator to first Dihedral in this Molecule.
      *
      * \param iterator on output, points to first Dihedral.
      */
      void begin(ConstDihedralIterator &iterator) const;
      #endif

      //@}

   private:

      /// Pointer to Species object.
      Species* speciesPtr_;

      /// Pointer to System containing this molecule, if any.
      /// This pointer should be null when the molecule is in no System.
      System* systemPtr_;

      /// Pointer to first atom in molecule.
      Atom* firstAtomPtr_;

      #ifdef SIMP_BOND
      /// Pointer to first bond in molecule.
      Bond* firstBondPtr_;
      #endif

      #ifdef SIMP_ANGLE
      /// Pointer to first angle in molecule.
      Angle* firstAnglePtr_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Pointer to first dihedral in molecule.
      Dihedral* firstDihedralPtr_;
      #endif

      /// Number of atoms in molecule.
      int nAtom_;

      #ifdef SIMP_BOND
      /// Number of bonds in molecule.
      int nBond_;
      #endif

      #ifdef SIMP_ANGLE
      /// Number of angles in molecule.
      int nAngle_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Number of dihedrals in molecule.
      int nDihedral_;
      #endif

      /// Integer index for this molecule within its Species
      int id_;

   // friends:

      class DeActivator;

   };

   // Inline member functions

   /*
   * Get the Species by reference.
   */
   inline Species& Molecule::species() const
   {
      assert(speciesPtr_);
      return *speciesPtr_;
   }

   /*
   * Get the parent System by reference.
   */
   inline System& Molecule::system() const
   {
      assert(systemPtr_);
      return *systemPtr_;
   }

   /*
   * Get number of atoms in this molecule.
   */
   inline int Molecule::nAtom() const
   {  return nAtom_; }

   #ifdef SIMP_BOND
   /*
   * Get number of bonds in this molecule.
   */
   inline int Molecule::nBond() const
   {  return nBond_; }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Get number of angles in this molecule.
   */
   inline int Molecule::nAngle() const
   {  return nAngle_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Get number of dihedrals in this molecule.
   */
   inline int Molecule::nDihedral() const
   {  return nDihedral_; }
   #endif

   /*
   * Get a specific Atom, referenced by an index.
   */
   inline Atom& Molecule::atom(int localIndex)
   {
      assert(firstAtomPtr_);
      assert(localIndex >= 0);
      assert(localIndex < nAtom_);
      return *(firstAtomPtr_ + localIndex);
   }

   /*
   * Get a specific Atom, referenced by an index.
   */
   inline const Atom& Molecule::atom(int localIndex) const
   {
      assert(firstAtomPtr_);
      assert(localIndex >= 0);
      assert(localIndex < nAtom_);
      return *(firstAtomPtr_ + localIndex);
   }

   #ifdef SIMP_BOND
   /*
   * Get a specific Bond by reference.
   */
   inline Bond& Molecule::bond(int localIndex)
   {
      assert(firstBondPtr_);
      assert(localIndex >= 0);
      assert(localIndex < nBond_);
      return *(firstBondPtr_ + localIndex);
   }

   /*
   * Get a specific Bond by const reference.
   */
   inline const Bond& Molecule::bond(int localIndex) const
   {
      assert(firstBondPtr_);
      assert(localIndex >= 0);
      assert(localIndex < nBond_);
      return *(firstBondPtr_ + localIndex);
   }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Get a specific Angle, referenced by an index.
   */
   inline Angle& Molecule::angle(int localIndex)
   {
      assert(firstAnglePtr_);
      assert(localIndex >= 0);
      assert(localIndex < nAngle_);
      return *(firstAnglePtr_ + localIndex);
   }

   /*
   * Get a specific Angle, referenced by an index.
   */
   inline const Angle& Molecule::angle(int localIndex) const
   {
      assert(firstAnglePtr_);
      assert(localIndex >= 0);
      assert(localIndex < nAngle_);
      return *(firstAnglePtr_ + localIndex);
   }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Get a specific Dihedral, referenced by an index.
   */
   inline Dihedral& Molecule::dihedral(int localIndex)
   {
      assert(firstDihedralPtr_);
      assert(localIndex >= 0);
      assert(localIndex < nDihedral_);
      return *(firstDihedralPtr_ + localIndex);
   }

   /*
   * Get a specific Dihedral, referenced by an index.
   */
   inline const Dihedral& Molecule::dihedral(int localIndex) const
   {
      assert(firstDihedralPtr_);
      assert(localIndex >= 0);
      assert(localIndex < nAngle_);
      return *(firstDihedralPtr_ + localIndex);
   }
   #endif

   /*
   * Get global index for this Molecule.
   */
   inline int Molecule::id() const
   {  return id_; }

   /*
   * Set AtomIterator to first Atom in this molecule.
   */
   inline void Molecule::begin(AtomIterator &iterator)
   {
      assert(firstAtomPtr_);
      assert(nAtom_ > 0);
      iterator.setCurrent(firstAtomPtr_);
      iterator.setEnd(firstAtomPtr_ + nAtom_);
   }

   /*
   * Set AtomIterator to first Atom in this molecule.
   */
   inline void Molecule::begin(ConstAtomIterator &iterator) const
   {
      assert(firstAtomPtr_);
      assert(nAtom_ > 0);
      iterator.setCurrent(firstAtomPtr_);
      iterator.setEnd(firstAtomPtr_ + nAtom_);
   }

   #ifdef SIMP_BOND
   /*
   * Set BondIterator to first Bond in this molecule.
   */
   inline void Molecule::begin(BondIterator &iterator)
   {
      assert(firstBondPtr_);
      assert(nBond_ > 0);
      iterator.setCurrent(firstBondPtr_);
      iterator.setEnd(firstBondPtr_ + nBond_);
   }

   /*
   * Set ConstBondIterator to first Bond in this molecule.
   */
   inline void Molecule::begin(ConstBondIterator &iterator) const
   {
      assert(firstBondPtr_);
      assert(nBond_ > 0);
      iterator.setCurrent(firstBondPtr_);
      iterator.setEnd(firstBondPtr_ + nBond_);
   }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Set AngleIterator to first Angle in this molecule.
   */
   inline void Molecule::begin(AngleIterator &iterator)
   {
      assert(firstAnglePtr_);
      assert(nAngle_ > 0);
      iterator.setCurrent(firstAnglePtr_);
      iterator.setEnd(firstAnglePtr_ + nAngle_);
   }

   /*
   * Set ConstAngleIterator to first Angle in this molecule.
   */
   inline void Molecule::begin(ConstAngleIterator &iterator) const
   {
      assert(firstAnglePtr_);
      assert(nAngle_ > 0);
      iterator.setCurrent(firstAnglePtr_);
      iterator.setEnd(firstAnglePtr_ + nAngle_);
   }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Set DihedralIterator to first Dihedral in this molecule.
   */
   inline void Molecule::begin(DihedralIterator &iterator)
   {
      assert(firstDihedralPtr_);
      assert(nDihedral_ > 0);
      iterator.setCurrent(firstDihedralPtr_);
      iterator.setEnd(firstDihedralPtr_ + nDihedral_);
   }

   /*
   * Set ConstDihedralIterator to first Dihedral in this molecule.
   */
   inline void Molecule::begin(ConstDihedralIterator &iterator) const
   {
      assert(firstDihedralPtr_);
      assert(nDihedral_ > 0);
      iterator.setCurrent(firstDihedralPtr_);
      iterator.setEnd(firstDihedralPtr_ + nDihedral_);
   }
   #endif

   /*
   * Set pointer to molecular Species.
   */
   inline void Molecule::setSpecies(Species &species)
   {  speciesPtr_ = &species; }

   /*
   * Set pointer to parent System.
   */
   inline void Molecule::setSystem(System &system)
   {  systemPtr_ = &system; }

   /*
   * Set the pointer to the parent System to null.
   */
   inline void Molecule::unsetSystem()
   {  systemPtr_ = 0; }

   /*
   * Set pointer to first Atom in molecule.
   */
   inline void Molecule::setFirstAtom(Atom &atom)
   {  firstAtomPtr_ = &atom; }

   /*
   * Set number of atoms per molecule.
   */
   inline void Molecule::setNAtom(int nAtom)
   {  nAtom_ = nAtom; }

   #ifdef SIMP_BOND
   /*
   * Set pointer to first Bond in molecule.
   */
   inline void Molecule::setFirstBond(Bond &bond)
   {  firstBondPtr_ = &bond; }

   /*
   * Set number of bonds per molecule.
   */
   inline void Molecule::setNBond(int nBond)
   {  nBond_ = nBond; }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Set pointer to first Angle in molecule.
   */
   inline void Molecule::setFirstAngle(Angle &angle)
   {  firstAnglePtr_ = &angle; }

   /*
   * Set number of angles per molecule.
   */
   inline void Molecule::setNAngle(int nAngle)
   {  nAngle_ = nAngle; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Set pointer to first Dihedral in molecule.
   */
   inline void Molecule::setFirstDihedral(Dihedral &dihedral)
   {  firstDihedralPtr_ = &dihedral; }

   /*
   * Set number of dihedrals per molecule.
   */
   inline void Molecule::setNDihedral(int nDihedral)
   {  nDihedral_ = nDihedral; }
   #endif

   /*
   * Set global index.
   */
   inline void Molecule::setId(int id)
   {  id_ = id; }

   // Inline member function of Atom class.

   /*
   * Return local id of an Atom within its parent molecule.
   *
   * This function is defined here, rather than in Atom.h, because its
   * implementation requires the inline method Molecule::atom(int).
   */
   inline int Atom::indexInMolecule() const
   {  return int(this - &molecule().atom(0));}

}
#endif
