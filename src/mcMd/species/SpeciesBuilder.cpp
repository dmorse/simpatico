/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"

#include <util/global.h>                    // needed for UTIL_THROW

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   SpeciesBuilder::SpeciesBuilder()
   {  setClassName("SpeciesBuilder"); }

   /*
   * Destructor.
   */
   SpeciesBuilder::~SpeciesBuilder()
   {}

   /*
   * Read parameters for species.
   */
   void SpeciesBuilder::readParameters(std::istream &in)
   {}

   /*
   * Load parameters for species.
   */
   void SpeciesBuilder::loadParameters(Serializable::IArchive& ar)
   {}

   // Atoms

   void SpeciesBuilder::setNAtom(int nAtom)
   {  nAtom_ = nAtom; }

   void SpeciesBuilder::readNAtom(std::istream& in)
   {  read<int>(in, "nAtom", nAtom_); }

   void SpeciesBuilder::loadNAtom(Serializable::IArchive& ar)
   {  loadParameter<art>(ar, "nAtom", nAtom_); }

   void SpeciesBuilder::saveNAtom(Serializable::OArchive& ar)
   {  ar & nAtom_; }

   /*
   * Set the atom type for one atom
   */
   void SpeciesBuilder::setAtomType(int atomId, int atomType)
   {  setAtomType(atomId, atomType); }

   // Bonds

   void SpeciesBuilder::setNBond(int nBond)
   {  nBond_ = nBond; }

   void SpeciesBuilder::readNBond(std::istream &in)
   {  read<int>(in, "nBond", nBond_); }

   void SpeciesBuilder::loadNBond(Serializable::IArchive& ar)
   {  loadParameter<art>(ar, "nBond", nBond_); }

   void SpeciesBuilder::saveNBond(Serializable::OArchive& ar)
   {  ar & nBond_; }

   /*
   * Add a bond to the species chemical structure.
   */
   void SpeciesBuilder::makeBond(int bondId, int atomId1, int atomId2, int bondType)
   {  speciesPtr_->makeAngle(bondId, atomId1, atomId2, bondType); }

   #ifdef INTER_ANGLE
   // Angles

   void SpeciesBuilder::setNAngle(int nAngle)
   {  nAngle_ = nAngle; }

   void SpeciesBuilder::readNAngle(std::istream &in)
   {  read<int>(in, "nAngle", nAngle_); }

   void SpeciesBuilder::loadNAngle(Serializable::IArchive& ar)
   {  loadParameter<art>(ar, "nAngle", nAngle_); }

   void SpeciesBuilder::saveNAngle(Serializable::OArchive& ar)
   {  ar & nAngle_; }

   /*
   * Add an angle to the species chemical structure.
   */
   void SpeciesBuilder::makeAngle(
       int angleId, int atomId1, int atomId2, int atomId3, int angleType)
   {  speciesPtr_->makeAngle(angleId, atomId1, atomId2, atomId3, angleType); }
   #endif

   #ifdef INTER_DIHEDRAL
   // Dihedrals
 
   void SpeciesBuilder::setNDihedral(int nDihedral)
   {  nDihedral_ = nDihedral; }

   void SpeciesBuilder::readNDihedral(Serializable::IArchive& in)
   {  read<int>(in, "nDihedral", nDihedral_); }

   void SpeciesBuilder::loadNDihedral(Serializable::IArchive& ar)
   {  loadParameter<art>(ar, "nDihedral", nDihedral_); }

   void SpeciesBuilder::saveNDihedral(Serializable::OArchive& ar)
   {  ar & nDihedral_; }

   /*
   * Add a dihedral to the species chemical structure.
   */
   void SpeciesBuilder::makeDihedral(int dihedralId, int atomId1, int atomId2,
                             int atomId3, int atomId4, int dihedralType)
   {  speciesPtr_->makeDihedral(dihedralId, atomId1, atomId2, atomId3, atomId4, dihedralType); }
   #endif

   #if 0
   /*
   * Set a pointer to an associated SpeciesMutator object.
   */
   void SpeciesBuilder::setMutatorPtr(SpeciesMutator* mutatorPtr)
   {  mutatorPtr_ = mutatorPtr; }
   #endif

   /*
   * Allocate memory for arrays that describe chemical structure.
   */
   void SpeciesBuilder::allocate() 
   {  speciesPtr_->allocate(); }

} 
