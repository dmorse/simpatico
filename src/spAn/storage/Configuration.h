#ifndef SPAN_STORAGE_H
#define SPAN_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>        // base class

#include <spAn/chemistry/Atom.h>              // member (template argument)
#include <spAn/chemistry/Group.h>             // member (template argument)
#include <spAn/chemistry/Species.h>           // member (template argument)
#include <spAn/storage/AtomStorage.h>         // member 
#include <spAn/storage/GroupStorage.h>        // member (template)

#include <util/boundary/Boundary.h>           // member 
#include <util/containers/DArray.h>           // member (template)
#include <util/containers/DSArray.h>          // member (template)
#include <util/containers/ArrayIterator.h>    // inline function


namespace SpAn 
{

   using namespace Util;

   /**
   * An instantaneous molecular dynamics configuration.
   *
   * A Configuration has:
   *   - a Util::Boundary
   *   - an AtomStorage container for atoms
   *   - a GroupStorage for each type of covalent group
   *
   * \ingroup SpAn_Storage_Module
   */
   class Configuration : public ParamComposite 
   {

   public:

      typedef ArrayIterator<Atom> AtomIterator;
      typedef ArrayIterator<Group <2> > BondIterator;
      typedef ArrayIterator<Group <3> > AngleIterator;
      typedef ArrayIterator<Group <4> > DihedralIterator;

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
      /// Get the Bond Configuration.
      GroupStorage<2>& bonds();
      #endif

      #ifdef INTER_ANGLE
      /// Get the Angle Configuration.
      GroupStorage<3>& angles();
      #endif

      #ifdef INTER_DIHEDRAL
      /// Get the Dihedral Configuration.
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
   inline Boundary& Configuration::boundary() 
   {  return boundary_; }

   inline AtomStorage& Configuration::atoms()
   {  return atoms_; }

   #ifdef INTER_BOND
   inline GroupStorage<2>& Configuration::bonds()
   {  return bonds_; }
   #endif

   #ifdef INTER_ANGLE
   inline GroupStorage<3>& Configuration::angles()
   {  return angles_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline GroupStorage<4>& Configuration::dihedrals()
   {  return dihedrals_; }
   #endif

   inline int Configuration::nSpecies() const
   {  return nSpecies_; }

   inline Species& Configuration::species(int i)
   {  return species_[i]; }

}
#endif
