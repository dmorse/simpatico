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
#include <mdPp/storage/GroupStorage.h>        // member (template)
#include <mdPp/chemistry/Species.h>           // member (template argument)

#include <util/boundary/Boundary.h>           // member 
#include <util/containers/DArray.h>           // member (template)
#include <util/containers/DSArray.h>          // member (template)
#include <util/containers/ArrayIterator.h>    // inline function


namespace MdPp 
{

   using namespace Util;

   /**
   * A post-processor for analyzing outputs of MD simulations.
   */
   class Storage : public ParamComposite 
   {

   public:

      typedef ArrayIterator<Atom> AtomIterator;
      typedef ArrayIterator<Group <2> > BondIterator;

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
  
      /**
      * Get the Boundary by non-const reference
      */
      Boundary& boundary();

      // Atom container interface

      /**
      * Return pointer to location for new atom.
      *
      * \param  global id for new atom
      * \return pointer to location of new atom
      */
      Atom* newAtomPtr();

      /**
      * Finalize addition of atom (allows lookup by id).
      */
      void addAtom();

      /**
      * Get a pointer to an atom by global id.
      */
      Atom* atomPtr(int id);

      /**
      * Initialize an iterator for atoms.
      */
      void initAtomIterator(AtomIterator& iter);

      /**
      * Get atom capacity (maximum id + 1).
      */ 
      int atomCapacity() const;

      /**
      * Get number of atoms.
      */ 
      int nAtom() const;

      // Group storage interface

      #ifdef INTER_BOND
      GroupStorage<2>& bonds();
      #endif

      #ifdef INTER_ANGLE
      GroupStorage<3>& angles();
      #endif

      #ifdef INTER_DIHEDRAL
      GroupStorage<4>& dihedrals();
      #endif

      // etc. for angles and dihedrals

      int nSpecies() const;

      Species& species(int i);

   private:
     
      /// Boundary object defines periodic boundary conditions.
      Boundary boundary_;

      /// Array of atom objects, added in order read from file.
      DSArray<Atom> atoms_;

      /// Pointers to atoms indexed by ids. Missing atoms are null pointers.
      DArray<Atom*> atomPtrs_;

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

   /*
   * Return number of atoms.
   */
   inline int Storage::nAtom() const
   {  return atoms_.size(); }

   /*
   * Get atom capacity (maximum id + 1).
   */ 
   inline
   int Storage::atomCapacity() const
   { return atoms_.capacity(); }

   /*
   * Return a pointer to an atom with a specific id.
   */
   inline Atom* Storage::atomPtr(int id)
   {  return atomPtrs_[id]; }

   /*
   * Initialize an iterator for atoms.
   */
   inline 
   void Storage::initAtomIterator(Storage::AtomIterator& iter)
   {  atoms_.begin(iter); }

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
