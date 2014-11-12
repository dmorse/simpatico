#ifndef MCMD_SPECIES_BUILDER_H
#define MCMD_SPECIES_BUILDER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <util/containers/ArrayStack.h>    // member template
#include <util/containers/DArray.h>        // member template
#include <util/containers/FSArray.h>       // member template

#include <mcMd/chemistry/SpeciesGroup.h>   // member template parameter

namespace McMd
{

   using namespace Util;

   class Species;

   /**
   * A SpeciesBuilder initializes a Species.
   *
   * \ingroup McMd_Species_Module
   */
   class Species : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      SpeciesBuilder();

      /**
      * Destructor.
      */
      virtual ~SpeciesBuilder();

      /**
      * Read parameters and initialize chemical structure for this species.
      *
      * \param in input stream.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

   protected:

      /**
      *
      * \param nBond number of bonds
      */
      void setNAtom(int nAtom);

      void readNAtom(std::istream& in);

      void loadNAtom(Serializable::IArchive& ar);

      void saveNAtom(Serializable::OArchive& ar);

      /**
      * Set the type for one atom in a generic molecule of this Species.
      *
      * \param atomId   local atom id, in range 0,..., nAtom - 1;
      * \param atomType atom type index
      */
      void setAtomType(int atomId, int atomType);

      // Bonds

      /**
      *
      *
      * \param nBond number of bonds
      */
      void setNBond(int nBond);

      /**
      *
      *
      * \param 
      */
      void readNBond(std::istream &in);

      /**
      *
      *
      * \param 
      */
      void loadNBond(Serializable::IArchive& ar);

      /**
      *
      *
      * \param 
      */
      void saveNBond(Serializable::OArchive& ar);

      /**
      * Add a bond to the chemical structure of a generic molecule.
      *
      * This function creates and adds a SpeciesBond object, and also adds
      * a reference to the list of bonds that are connected to each atom.
      *
      * \param bondId   local index of bond within a molecule
      * \param atomId1  local index of 1st atom in bond
      * \param atomId2  local index of 2nd atom in bond
      * \param bondType bond type index
      */
      void makeBond(int bondId, int atomId1, int atomId2, int bondType);

      #ifdef INTER_ANGLE
      /**
      *
      * \param nAngle number of angles
      */
      void setNAngle(int nAngle);

      /**
      *
      *
      * \param 
      */
      void readNAngle(std::istream &in);

      /**
      *
      *
      * \param 
      */
      void loadNAngle(Serializable::IArchive& ar);

      /**
      *
      *
      * \param 
      */
      void saveNAngle(Serializable::OArchive& ar);

      /**
      * Add an angle to the chemical structure of a generic molecule.
      *
      * This function creates and adds a SpeciesAngle object, and also adds
      * a reference to the list of angles that are connected to each atom.
      *
      * \param angleId   local index of angle within a molecule
      * \param atomId1   local index of 1st atom in angle
      * \param atomId2   local index of 2nd atom in angle
      * \param atomId3   local index of 3rd atom in angle
      * \param angleType angle type index
      */
      void makeAngle(int angleId, int atomId1, int atomId2, int atomId3,
                     int angleType);
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      *
      * \param nDihedral number of dihedrals
      */
      void setNDihedral(int nDihedral);

      /**
      *
      *
      * \param 
      */
      void readNDihedral(Serializable::IArchive& in);

      /**
      *
      *
      * \param 
      */
      void loadNDihedral(Serializable::IArchive& ar);

      /**
      *
      *
      * \param 
      */
      void saveNDihedral(Serializable::OArchive& ar);

      /**
      * Add a dihedral to the chemical structure of a generic molecule.
      *
      * This function creates and adds a SpeciesDihedral object, and also adds
      * a reference to the list of dihedrals that are connected to each atom.
      *
      * \param dihedralId  local index of dihedral within a molecule
      * \param atomId1     local index of 1st atom in dihedral
      * \param atomId2     local index of 2nd atom in dihedral
      * \param atomId3     local index of 3rd atom in dihedral
      * \param atomId4     local index of 4th atom in dihedral
      * \param dihedralType dihedral type index
      */
      void makeDihedral(int dihedralId, int atomId1, int atomId2, int atomId3,
                     int atomId4, int dihedralType);
      #endif

      /**
      * Allocate chemical structure arrays.
      *
      * This function allocates the arrays that are used to define the
      * chemical structure of a generic molecule, such as atomTypeIds_,
      * speciesBonds_, atomBondIdArrays_, speciesAngles_, etc.
      *
      * Precondition: nAtom_, nBond_, nAngles_, etc. must be set before
      * entry.
      */
      void allocate();

   private:

      Species* speciesPtr_;

   }

}
#endif
