#ifndef MCMD_SYSTEM_INTERFACE_H
#define MCMD_SYSTEM_INTERFACE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>               // parent (provides typedefs)
#include <util/boundary/Boundary.h>               // member
#include <mcMd/chemistry/Molecule.h>              // member template parameter

#include <util/containers/DArray.h>               // member template
#include <util/containers/PArrayIterator.h>       // inline function begin()
#include <util/global.h>

#include <iostream>
#include <string>

namespace McMd
{

   class Simulation;

   using namespace Util;

   /**
   * An interface to a System.
   *
   * A SystemInterface is an interface for a system that is intended
   * to be used as a protected or private base class for classes 
   * that evaluate properties of a system. It provides functions
   * to directly access:
   *
   *  - an ArraySet of Molecule objects of each Species.
   *  - a Boundary.
   *
   * \ingroup McMd_System_Module
   */
   class SystemInterface
   {

   public:

      /**
      * Constructor.
      *
      * \param parent parent System
      */
      SystemInterface(System& parent);

      /**
      * Destructor.
      */
      virtual ~SystemInterface();

      /**
      * Get the parent System by reference.
      */
      System& system() const;

   protected:

      /// \name Accessors (Miscellaneous)
      //@{

      /// Get the Boundary by reference.
      Boundary& boundary() const;

      /// Get the parent Simulation by reference.
      Simulation& simulation() const;

      //@}
      /// \name Molecule Set Accessors
      //@{

      /**
      * Get the number of molecules of one Species in this SystemInterface.
      *
      * \param speciesId integer Id for a Species.
      * \return number of molecules of specified Species in this SystemInterface.
      */
      int nMolecule(int speciesId) const;

      /**
      * Return the total number of atoms in this SystemInterface.
      */
      int nAtom() const;

      /**
      * Is this an empty SystemInterface (i.e., one with no molecules) ?
      */
      bool isEmpty() const;

      /**
      * Initialize an iterator for molecules of one species in this SystemInterface.
      *
      * \param speciesId integer Id for the desired Species (input)
      * \param iterator  molecule iterator (output)
      */
      void begin(int speciesId, System::MoleculeIterator& iterator);

      /**
      * Initialize a const iterator for molecules of one species in this SystemInterface.
      *
      * \param speciesId integer Id for the desired Species (input)
      * \param iterator  molecule iterator (output)
      */
      void begin(int speciesId, System::ConstMoleculeIterator& iterator) const;

      //@}
      /// \name Potential Energy Queries
      //@{

      #ifdef SIMP_BOND
      /// Does a bond potential exist?
      bool hasBonds() const;
      #endif

      #ifdef SIMP_ANGLE
      /// Does an angle potential exist?
      bool hasAngles() const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Does a dihedral potential exist?
      bool hasDihedrals() const;
      #endif

      #ifdef MCMD_LINK
      /// Does a link potential exist?
      bool hasLinks() const;
      #endif

      #ifdef SIMP_EXTERNAL
      /// Does an external potential exist?
      bool hasExternal() const;
      #endif

      #ifdef SIMP_TETHER
      /// Does a tether potential exist?
      bool hasTethers() const;
      #endif

      //@}

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

      /// Pointer to parent System.
      System* systemPtr_;

      /**
      * Pointer to DArray containing one System::MoleculeSet for each Species.
      *
      * The MoleculeSet (*moleculeSetsPtr_)[i] contains the molecules in
      * this SystemInterface that belong to Species i of the parent simulation.
      */
      DArray< System::MoleculeSet >* moleculeSetsPtr_;

      /// Pointer to Boundary object for actual boundary.
      Boundary*   boundaryPtr_;

      #ifdef SIMP_BOND
      // Does a bond potential exist?
      bool hasBonds_;
      #endif

      #ifdef SIMP_ANGLE
      // Does an angle potential exist?
      bool hasAngles_;
      #endif

      #ifdef SIMP_DIHEDRAL
      // Does a dihedral potential exist?
      bool hasDihedrals_;
      #endif

      #ifdef MCMD_LINK
      // Does a link potential exist?
      bool hasLinks_;
      #endif

      #ifdef SIMP_EXTERNAL
      // Does an external potential exist?
      bool hasExternal_;
      #endif

      #ifdef SIMP_TETHER
      // Does a tether potential exist?
      bool hasTethers_;
      #endif

   };

   // Inline functions

   /*
   * Get the parent Simulation by reference.
   */
   inline Simulation& SystemInterface::simulation() const
   {
      assert(simulationPtr_);
      return *simulationPtr_;
   }

   /*
   * Get the parent System by reference.
   */
   inline System& SystemInterface::system() const
   {
      assert(systemPtr_);
      return *systemPtr_;
   }

   /*
   * Get the Boundary by reference.
   */
   inline Boundary& SystemInterface::boundary() const
   {
      assert(boundaryPtr_);
      return *boundaryPtr_;
   }

   /*
   * Get the number of molecules of a specific Species in this SystemInterface.
   */
   inline int SystemInterface::nMolecule(int speciesId) const
   {
      assert(moleculeSetsPtr_);
      return (*moleculeSetsPtr_)[speciesId].size();
   }

   /*
   * Initialize a System::MoleculeIterator for molecules of one Species.
   */
   inline void
   SystemInterface::begin(int speciesId, System::MoleculeIterator& iterator)
   {
      assert(moleculeSetsPtr_);
      (*moleculeSetsPtr_)[speciesId].begin(iterator);
   }

   /*
   * Initialize a System::ConstMoleculeIterator for molecules of one Species.
   */
   inline void
   SystemInterface::begin(int speciesId, System::ConstMoleculeIterator& iterator) const
   {
      assert(moleculeSetsPtr_);
      (*moleculeSetsPtr_)[speciesId].begin(iterator);
   }

   #ifdef SIMP_BOND
   /*
   * Does a bond potential exist?
   */
   inline bool SystemInterface::hasBonds() const
   {  return hasBonds_; }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Does an angle potential exist?
   */
   inline bool SystemInterface::hasAngles() const
   {  return hasAngles_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /// Does a dihedral potential exist?
   inline bool SystemInterface::hasDihedrals() const
   {  return hasDihedrals_; }
   #endif

   #ifdef MCMD_LINK
   /// Does a link potential exist?
   inline bool SystemInterface::hasLinks() const
   { return hasLinks_; }
   #endif

   #ifdef SIMP_EXTERNAL
   /// Does an external potential exist?
   inline bool SystemInterface::hasExternal() const
   { return hasExternal_; }
   #endif

   #ifdef SIMP_TETHER
   /// Does a tether potential exist?
   inline bool SystemInterface::hasTethers() const
   { return hasTethers_; }
   #endif

}
#endif
