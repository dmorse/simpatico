#ifndef DDMD_SP_STORAGE_H
#define DDMD_SP_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>        // base class

#include <ddMd/sp/chemistry/SpAtom.h>              // member (template argument)
#include <ddMd/sp/chemistry/SpGroup.h>             // member (template argument)
#include <ddMd/sp/chemistry/SpSpecies.h>           // member (template argument)
#include <ddMd/sp/storage/SpAtomStorage.h>         // member 
#include <ddMd/sp/storage/SpGroupStorage.h>        // member (template)

#include <util/boundary/Boundary.h>           // member 
#include <util/containers/DArray.h>           // member (template)
#include <util/containers/DSArray.h>          // member (template)
#include <util/containers/ArrayIterator.h>    // inline function


namespace DdMd 
{

   using namespace Util;

   /**
   * A snapshot of a molecular dynamics configuration configuration.
   *
   * A SpConfiguration has:
   *   - a Boundary
   *   - an SpAtomStorage container for atoms
   *   - a SpGroupStorage for each type of covalent group
   *
   * \ingroup DdMd_Sp_Storage_Module
   */
   class SpConfiguration : public ParamComposite 
   {

   public:

      typedef ArrayIterator<SpAtom> AtomIterator;
      typedef ArrayIterator<SpGroup <2> > BondIterator;
      typedef ArrayIterator<SpGroup <3> > AngleIterator;
      typedef ArrayIterator<SpGroup <4> > DihedralIterator;

      /**
      * Constructor
      */
      SpConfiguration();

      /**
      * Destructor
      */
      ~SpConfiguration();

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

      /// Get the SpAtomStorage.
      SpAtomStorage& atoms();

      #ifdef INTER_BOND
      /// Get the Bond SpConfiguration.
      SpGroupStorage<2>& bonds();
      #endif

      #ifdef INTER_ANGLE
      /// Get the Angle SpConfiguration.
      SpGroupStorage<3>& angles();
      #endif

      #ifdef INTER_DIHEDRAL
      /// Get the Dihedral SpConfiguration.
      SpGroupStorage<4>& dihedrals();
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
      SpSpecies& species(int i);

   private:
     
      /// Boundary object defines periodic boundary conditions.
      Boundary boundary_;

      /// SpAtomStorage object.
      SpAtomStorage atoms_;

      #ifdef INTER_BOND
      /// Array of bond objects, added in order read from file.
      SpGroupStorage<2> bonds_;
      #endif

      #ifdef INTER_ANGLE
      /// Array of angle objects, added in order read from file.
      SpGroupStorage<3> angles_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Array of dihedral objects, added in order read from file.
      SpGroupStorage<4> dihedrals_;
      #endif

      /// Array of SpSpecies objects.
      DArray<SpSpecies> species_;

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
   inline Boundary& SpConfiguration::boundary() 
   {  return boundary_; }

   inline SpAtomStorage& SpConfiguration::atoms()
   {  return atoms_; }

   #ifdef INTER_BOND
   inline SpGroupStorage<2>& SpConfiguration::bonds()
   {  return bonds_; }
   #endif

   #ifdef INTER_ANGLE
   inline SpGroupStorage<3>& SpConfiguration::angles()
   {  return angles_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline SpGroupStorage<4>& SpConfiguration::dihedrals()
   {  return dihedrals_; }
   #endif

   inline int SpConfiguration::nSpecies() const
   {  return nSpecies_; }

   inline SpSpecies& SpConfiguration::species(int i)
   {  return species_[i]; }

}
#endif
