#ifndef MCMD_SUB_SYSTEM_H
#define MCMD_SUB_SYSTEM_H

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
   * A SubSystem is a partial copy of a System that has:
   *
   *  - an ArraySet of Molecule objects of each Species.
   *  - a Boundary.
   *
   * A SubSystem must be associated with a parent System.
   *
   * SubSystem is designed as a base class for classes that
   * evaluate configurational properties, energy, and forces.
   *
   * \ingroup McMd_System_Module
   */
   class SubSystem
   {

   public:

      /**
      * Constructor.
      *
      * \param parent parent System
      */
      SubSystem(System& parent);

      /**
      * Destructor.
      */
      virtual ~SubSystem();

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
      * Get the number of molecules of one Species in this SubSystem.
      *
      * \param speciesId integer Id for a Species.
      * \return number of molecules of specified Species in this SubSystem.
      */
      int nMolecule(int speciesId) const;

      /**
      * Return the total number of atoms in this SubSystem.
      */
      int nAtom() const;

      /**
      * Is this an empty SubSystem (i.e., one with no molecules) ?
      */
      bool isEmpty() const;

      /**
      * Initialize an iterator for molecules of one species in this SubSystem.
      *
      * \param speciesId integer Id for the desired Species (input)
      * \param iterator  molecule iterator (output)
      */
      void begin(int speciesId, System::MoleculeIterator& iterator);

      /**
      * Initialize a const iterator for molecules of one species in this SubSystem.
      *
      * \param speciesId integer Id for the desired Species (input)
      * \param iterator  molecule iterator (output)
      */
      void begin(int speciesId, System::ConstMoleculeIterator& iterator) const;

      //@}
      /// \name Potential Energy Queries
      //@{

      #ifdef INTER_BOND
      /// Does a bond potential exist?
      bool hasBonds() const;
      #endif

      #ifdef INTER_ANGLE
      /// Does an angle potential exist?
      bool hasAngles() const;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Does a dihedral potential exist?
      bool hasDihedrals() const;
      #endif

      #ifdef MCMD_LINK
      /// Does a link potential exist?
      bool hasLinks() const;
      #endif

      #ifdef INTER_EXTERNAL
      /// Does an external potential exist?
      bool hasExternal() const;
      #endif

      #ifdef INTER_TETHER
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
      * this SubSystem that belong to Species i of the parent simulation.
      */
      DArray< System::MoleculeSet >* moleculeSetsPtr_;

      /// Pointer to Boundary object for actual boundary.
      Boundary*   boundaryPtr_;

      #ifdef INTER_BOND
      // Does a bond potential exist?
      bool hasBonds_;
      #endif

      #ifdef INTER_ANGLE
      // Does an angle potential exist?
      bool hasAngles_;
      #endif

      #ifdef INTER_DIHEDRAL
      // Does a dihedral potential exist?
      bool hasDihedrals_;
      #endif

      #ifdef MCMD_LINK
      // Does a link potential exist?
      bool hasLinks_;
      #endif

      #ifdef INTER_EXTERNAL
      // Does an external potential exist?
      bool hasExternal_;
      #endif

      #ifdef INTER_TETHER
      // Does a tether potential exist?
      bool hasTethers_;
      #endif

   };

   // Inline functions

   /*
   * Get the parent Simulation by reference.
   */
   inline Simulation& SubSystem::simulation() const
   {
      assert(simulationPtr_);
      return *simulationPtr_;
   }

   /*
   * Get the parent System by reference.
   */
   inline System& SubSystem::system() const
   {
      assert(systemPtr_);
      return *systemPtr_;
   }

   /*
   * Get the Boundary by reference.
   */
   inline Boundary& SubSystem::boundary() const
   {
      assert(boundaryPtr_);
      return *boundaryPtr_;
   }

   /*
   * Get the number of molecules of a specific Species in this SubSystem.
   */
   inline int SubSystem::nMolecule(int speciesId) const
   {
      assert(moleculeSetsPtr_);
      return (*moleculeSetsPtr_)[speciesId].size();
   }

   /*
   * Initialize a System::MoleculeIterator for molecules of one Species.
   */
   inline void
   SubSystem::begin(int speciesId, System::MoleculeIterator& iterator)
   {
      assert(moleculeSetsPtr_);
      (*moleculeSetsPtr_)[speciesId].begin(iterator);
   }

   /*
   * Initialize a System::ConstMoleculeIterator for molecules of one Species.
   */
   inline void
   SubSystem::begin(int speciesId, System::ConstMoleculeIterator& iterator) const
   {
      assert(moleculeSetsPtr_);
      (*moleculeSetsPtr_)[speciesId].begin(iterator);
   }

   #ifdef INTER_BOND
   /*
   * Does a bond potential exist?
   */
   inline bool SubSystem::hasBonds() const
   {  return hasBonds_; }
   #endif

   #ifdef INTER_ANGLE
   /*
   * Does an angle potential exist?
   */
   inline bool SubSystem::hasAngles() const
   {  return hasAngles_; }
   #endif

   #ifdef INTER_DIHEDRAL
   /// Does a dihedral potential exist?
   inline bool SubSystem::hasDihedrals() const
   {  return hasDihedrals_; }
   #endif

   #ifdef MCMD_LINK
   /// Does a link potential exist?
   inline bool SubSystem::hasLinks() const
   { return hasLinks_; }
   #endif

   #ifdef INTER_EXTERNAL
   /// Does an external potential exist?
   inline bool SubSystem::hasExternal() const
   { return hasExternal_; }
   #endif

   #ifdef INTER_TETHER
   /// Does a tether potential exist?
   inline bool SubSystem::hasTethers() const
   { return hasTethers_; }
   #endif

}
#endif
