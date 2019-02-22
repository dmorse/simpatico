#ifndef MDPP_CONFIGURATION_H
#define MDPP_CONFIGURATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class

#include <mdPp/chemistry/Atom.h>             // member (template argument)
#include <mdPp/chemistry/Group.h>            // member (template argument)
#include <mdPp/storage/AtomStorage.h>        // member 
#include <mdPp/storage/GroupStorage.h>       // member (template)
#include <mdPp/storage/SpeciesStorage.h>     // member (template argument)
#include <simp/boundary/Boundary.h>          // member 
#include <util/containers/DArray.h>          // member (template)
#include <util/containers/DSArray.h>         // member (template)
#include <util/containers/ArrayIterator.h>   // inline function

namespace MdPp 
{

   using namespace Util;
   using namespace Simp;

   /**
   * An instantaneous molecular dynamics configuration.
   *
   * A Configuration has:
   *   - a Simp::Boundary
   *   - an AtomStorage container for atoms
   *   - a GroupStorage for each type of covalent group
   *   - a SpeciesStorage for each species (optionally)
   *
   * \ingroup MdPp_Storage_Module
   */
   class Configuration : public ParamComposite 
   {

   public:

      typedef ArrayIterator<Atom> AtomIterator;
      typedef ArrayIterator<Group <2> > BondIterator;
      typedef ArrayIterator<Group <3> > AngleIterator;
      typedef ArrayIterator<Group <4> > DihedralIterator;
      typedef ArrayIterator<Group <4> > ImproperIterator;

      /**
      * Constructor
      */
      Configuration();

      /**
      * Destructor
      */
      ~Configuration();

      using ParamComposite::readParam;

      /**
      * Open, read, and close parameter file.
      */
      void readParam(const char* filename);

      /**
      * Read parameters.
      */
      void readParameters(std::istream& in);

      // Accessors for members (non-const reference)

      /**
      * Get the Boundary by reference
      */
      Boundary& boundary();

      /**
      * Number of species.
      *
      * If nSpecies == 0, all species and molecule info is disabled.
      */
      int nSpecies() const;

      /**
      * Get a particular SpeciesStorage identified by index.
      *
      * \throw Exception if nSpecies == 0 or i is out of bounds.
      *
      * \param i species index
      */
      SpeciesStorage& species(int i);

      /**
      * Get the AtomStorage by reference.
      */
      AtomStorage& atoms();

      #ifdef SIMP_BOND
      /**
      * Get the Bond storage by reference.
      */
      GroupStorage<2>& bonds();
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get the Angle storage by reference.
      */
      GroupStorage<3>& angles();
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Get the Dihedral storage by reference.
      */
      GroupStorage<4>& dihedrals();

      /**
      * Get the Improper dihedral storage by reference.
      */
      GroupStorage<4>& impropers();
      #endif

      // Initialization

      /**
      * Allocate space for species information.
      *
      * Call in configuration file reader before reading species
      * structure information.
      *
      * \param nSpecies number of species in system
      */
      void setNSpecies(int nSpecies);

      /**
      * Clear all atoms and groups.
      */
      void clear();
  
      /**
      * Set atom context data for all atoms, assuming consecutive ids.
      * 
      * \return true for normal completion, false if error is detected.
      */
      void setAtomContexts();

      /**
      * Add all atoms in AtomStorage to the SpeciesStorage.
      */
      void addAtomsToSpecies();

      /*
      * Create all covalent groups from species templates.
      */
      void makeGroups();

   private:
     
      /// Boundary object defines periodic boundary conditions.
      Boundary boundary_;

      /// Array of Species objects.
      DArray<SpeciesStorage> species_;

      /// AtomStorage object.
      AtomStorage atoms_;

      #ifdef SIMP_BOND
      /// Array of bond objects, added in order read from file.
      GroupStorage<2> bonds_;
      #endif

      #ifdef SIMP_ANGLE
      /// Array of angle objects, added in order read from file.
      GroupStorage<3> angles_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Array of dihedral objects, added in order read from file.
      GroupStorage<4> dihedrals_;

      /// Array of improper objects, added in order read from file.
      GroupStorage<4> impropers_;
      #endif

      /// Number of species (set to zero to disable)
      int nSpecies_;

      // Variables atomCapacity, bondCapacity etc. are only used to store
      // values read from a parameter file. Should not be accessed in any
      // other way.

      /// Maximum number of atoms = max id + 1 (used to allocate arrays).
      int atomCapacity_;

      #ifdef SIMP_BOND
      /// Maximum number of bonds (used to allocate array).
      int bondCapacity_;
      #endif

      #ifdef SIMP_ANGLE
      /// Maximum number of angles (used to allocate array).
      int angleCapacity_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Maximum number of dihedrals (used to allocate array).
      int dihedralCapacity_;
      #endif

      #ifdef SIMP_IMPROPER
      /// Maximum number of impropers (used to allocate array).
      int improperCapacity_;
      #endif

      #ifdef SIMP_BOND
      /*
      * Create all bonds from species templates.
      */
      void makeBonds();
      #endif

      #ifdef SIMP_ANGLE
      /*
      * Create all angles from species templates.
      */
      void makeAngles();
      #endif

      #ifdef SIMP_DIHEDRAL
      /*
      * Create all dihedrals from species templates.
      */
      void makeDihedrals();
      #endif

      /*
      * Create all Group<N> objects for one species.
      */
      template <int N>
      void makeSpeciesGroups(
          GroupStorage<N>& storage,
          const DArray< SpeciesGroup<N> >& speciesGroups,
          int nMol, int nAtom, int nGroup, 
          int& firstAtomId, int& groupId);

   };

   // Inline functions

   /*
   * Return the Boundary by reference.
   */
   inline Boundary& Configuration::boundary() 
   {  return boundary_; }

   inline AtomStorage& Configuration::atoms()
   {  return atoms_; }

   #ifdef SIMP_BOND
   inline GroupStorage<2>& Configuration::bonds()
   {  return bonds_; }
   #endif

   #ifdef SIMP_ANGLE
   inline GroupStorage<3>& Configuration::angles()
   {  return angles_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   inline GroupStorage<4>& Configuration::dihedrals()
   {  return dihedrals_; }

   inline GroupStorage<4>& Configuration::impropers()
   {  return impropers_; }
   #endif

   inline int Configuration::nSpecies() const
   {  return nSpecies_; }

   inline SpeciesStorage& Configuration::species(int i)
   {
      UTIL_CHECK(species_.isAllocated());
      return species_[i]; 
   }

}
#endif
