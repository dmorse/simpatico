#ifndef MDPP_STORAGE_H
#define MDPP_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>        // base class

#include <mdPp/chemistry/Atom.h>              // member (template argument)
#include <mdPp/chemistry/Group.h>             // member (template argument)
#include <mdPp/chemistry/Species.h>           // member (template argument)
#include <mdPp/storage/AtomStorage.h>         // member 
#include <mdPp/storage/GroupStorage.h>        // member (template)

#include <util/boundary/Boundary.h>           // member 
#include <util/containers/DArray.h>           // member (template)
#include <util/containers/DSArray.h>          // member (template)
#include <util/containers/ArrayIterator.h>    // inline function


namespace MdPp 
{

   using namespace Util;

   /**
   * A molecular dynamics system configuration.
   *
   * A Storage has:
   *   - a Boundary
   *   - an array of atoms
   *   - a GroupStorage for each type of covalent group
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

      AtomStorage& atoms();

      #ifdef INTER_BOND
      GroupStorage<2>& bonds();
      #endif

      #ifdef INTER_ANGLE
      GroupStorage<3>& angles();
      #endif

      #ifdef INTER_DIHEDRAL
      GroupStorage<4>& dihedrals();
      #endif

      int nSpecies() const;

      Species& species(int i);

   private:
     
      /// Boundary object defines periodic boundary conditions.
      Boundary boundary_;

      /// AtomStorage object
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

      /// Pointer to new atom.
      Atom* newAtomPtr_;

      /// Maximum allowed atom id + 1 (used to allocate arrays).
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

   // inline functions

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
