#ifndef MDCF_STORAGE_H
#define MDCF_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>        // base class

#include <mdCf/chemistry/Atom.h>              // member (template argument)
#include <mdCf/chemistry/Group.h>             // member (template argument)
#include <mdCf/chemistry/Species.h>           // member (template argument)
#include <mdCf/storage/AtomStorage.h>         // member 
#include <mdCf/storage/GroupStorage.h>        // member (template)

#include <util/boundary/Boundary.h>           // member 
#include <util/containers/DArray.h>           // member (template)
#include <util/containers/DSArray.h>          // member (template)
#include <util/containers/ArrayIterator.h>    // inline function


namespace MdCf 
{

   using namespace Util;

   /**
   * A snapshot of a molecular dynamics system configuration.
   *
   * A Storage has:
   *   - a Boundary
   *   - an AtomStorage container for atoms
   *   - a GroupStorage for each type of covalent group
   *
   * \ingroup MdCf_Storage_Module
   */
   class Storage : public ParamComposite 
   {

   public:

      typedef ArrayIterator<Atom> AtomIterator;
      typedef ArrayIterator<Group <2> > BondIterator;
      typedef ArrayIterator<Group <3> > AngleIterator;
      typedef ArrayIterator<Group <4> > DihedralIterator;

      /**
      * Constructor
      */
      Storage();

      /**
      * Destructor
      */
      ~Storage();

      using ParamComposite::readParam;

      /**
      * Open, read, and close parameter file.
      */
      void readParam(const char* filename);

      /**
      * Read parameters.
      */
      void readParameters(std::istream& in);

      /**
      * Clear all atoms and groups.
      */
      void clear();
  
      // Accessors for members (non-const reference)

      /**
      * Get the Boundary by non-const reference
      */
      Boundary& boundary();

      /// Get the AtomStorage.
      AtomStorage& atoms();

      #ifdef INTER_BOND
      /// Get the Bond Storage.
      GroupStorage<2>& bonds();
      #endif

      #ifdef INTER_ANGLE
      /// Get the Angle Storage.
      GroupStorage<3>& angles();
      #endif

      #ifdef INTER_DIHEDRAL
      /// Get the Dihedral Storage.
      GroupStorage<4>& dihedrals();
      #endif

      /**
      * Number of species.
      *
      * If nSpecies == 0, all species and molecule info is disabled.
      */
      int nSpecies() const;

      /**
      * Get a particular species identified by index.
      *
      * \param i species index
      */
      Species& species(int i);

   private:
     
      /// Boundary object defines periodic boundary conditions.
      Boundary boundary_;

      /// AtomStorage object.
      AtomStorage atoms_;

      #ifdef INTER_BOND
      /// Array of bond objects, added in order read from file.
      GroupStorage<2> bonds_;
      #endif

      #ifdef INTER_ANGLE
      /// Array of angle objects, added in order read from file.
      GroupStorage<3> angles_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Array of dihedral objects, added in order read from file.
      GroupStorage<4> dihedrals_;
      #endif

      /// Array of Species objects.
      DArray<Species> species_;

      /// Maximum number of atoms = max id + 1 (used to allocate arrays).
      int atomCapacity_;

      #ifdef INTER_BOND
      /// Maximum number of bonds (used to allocate array).
      int bondCapacity_;
      #endif

      #ifdef INTER_ANGLE
      /// Maximum number of angles (used to allocate array).
      int angleCapacity_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Maximum number of dihedrals (used to allocate array).
      int dihedralCapacity_;
      #endif

      /// Number of species (set to zero to disable)
      int nSpecies_;

   };

   // Inline functions

   /*
   * Return the Boundary by reference.
   */
   inline Boundary& Storage::boundary() 
   {  return boundary_; }

   inline AtomStorage& Storage::atoms()
   {  return atoms_; }

   #ifdef INTER_BOND
   inline GroupStorage<2>& Storage::bonds()
   {  return bonds_; }
   #endif

   #ifdef INTER_ANGLE
   inline GroupStorage<3>& Storage::angles()
   {  return angles_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline GroupStorage<4>& Storage::dihedrals()
   {  return dihedrals_; }
   #endif

   inline int Storage::nSpecies() const
   {  return nSpecies_; }

   inline Species& Storage::species(int i)
   {  return species_[i]; }

}
#endif
