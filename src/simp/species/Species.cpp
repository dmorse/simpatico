/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"

#include <util/global.h>                    
#include <simp/species/SpeciesGroup.tpp>  

namespace Simp
{

   using namespace Util;

   /*
   * Default constructor.
   */
   Species::Species()
    : id_(-1),
      moleculeCapacity_(0),
      nAtom_(0),
      atomTypeIds_(),
      #ifdef SIMP_BOND
      nBond_(0),
      speciesBonds_(),
      #endif
      #ifdef SIMP_ANGLE
      nAngle_(0),
      speciesAngles_(),
      #endif
      #ifdef SIMP_DIHEDRAL
      nDihedral_(0),
      speciesDihedrals_(),
      #endif
      mutatorPtr_(0)
   {  setClassName("Species"); }

   /*
   * Destructor.
   */
   Species::~Species()
   {}

   /*
   * Read parameters for species.
   */
   void Species::readParameters(std::istream &in)
   {
      read<int>(in, "moleculeCapacity", moleculeCapacity_);

      // Define chemical structure, after reading any parameters.
      readSpeciesParam(in);

      // Check validity, throw exception if not valid
      isValid();
   }

   /*
   * Load parameters for species.
   */
   void Species::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<int>(ar, "moleculeCapacity", moleculeCapacity_);

      // Read and define chemical structure.
      loadSpeciesParam(ar);

      // Check validity, throw exception if not valid
      isValid();
   }

   /*
   * Define molecular structure (default implementation).
   */
   void Species::readSpeciesParam(std::istream &in)
   {
      read<int>(in, "nAtom", nAtom_);
      #ifdef SIMP_BOND
      read<int>(in, "nBond", nBond_);
      #endif
      #ifdef SIMP_ANGLE
      read<int>(in, "nAngle", nAngle_);
      #endif
      #ifdef SIMP_DIHEDRAL
      read<int>(in, "nDihedral", nDihedral_);
      #endif
      allocate();

      readDArray<int>(in, "atomTypeIds", atomTypeIds_, nAtom_);

      #ifdef SIMP_BOND
      if (nBond_ > 0) {
         readDArray<SpeciesBond>(in, "speciesBonds", speciesBonds_, 
                                 nBond_);
         // Make atomBondIdArrays
         int bondId, atomId1, atomId2;
         for (bondId = 0; bondId < nBond_; ++bondId) {
            atomId1 = speciesBond(bondId).atomId(0);
            atomId2 = speciesBond(bondId).atomId(1);
            atomBondIdArrays_[atomId1].append(bondId);
            atomBondIdArrays_[atomId2].append(bondId);
         }
      }
      #endif

      #ifdef SIMP_ANGLE
      if (nAngle_ > 0) {
         readDArray<SpeciesAngle>(in, "speciesAngles", speciesAngles_,
                                  nAngle_);
         // Make atomAngleIdArrays
         int angleId, atomId1, atomId2, atomId3;
         for (angleId = 0; angleId < nAngle_; ++angleId) {
            atomId1 = speciesAngle(angleId).atomId(0);
            atomId2 = speciesAngle(angleId).atomId(1);
            atomId3 = speciesAngle(angleId).atomId(2);
            atomAngleIdArrays_[atomId1].append(angleId);
            atomAngleIdArrays_[atomId2].append(angleId);
            atomAngleIdArrays_[atomId3].append(angleId);
         }
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      // Make atomDihedralIdArrays
      if (nDihedral_ > 0) {
         readDArray<SpeciesDihedral>(in, "speciesDihedrals", speciesDihedrals_,
                                     nDihedral_);
         int dihedralId;
         int tAtomId1, tAtomId2, tAtomId3, tAtomId4;
         for (dihedralId = 0; dihedralId < nDihedral_; ++dihedralId) {
            tAtomId1 = speciesDihedral(dihedralId).atomId(0);
            tAtomId2 = speciesDihedral(dihedralId).atomId(1);
            tAtomId3 = speciesDihedral(dihedralId).atomId(2);
            tAtomId4 = speciesDihedral(dihedralId).atomId(3);
            atomDihedralIdArrays_[tAtomId1].append(dihedralId);
            atomDihedralIdArrays_[tAtomId2].append(dihedralId);
            atomDihedralIdArrays_[tAtomId3].append(dihedralId);
            atomDihedralIdArrays_[tAtomId4].append(dihedralId);
         }
      }
      #endif

   }

   /*
   * Load molecular structure (default implementation).
   */
   void Species::loadSpeciesParam(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar, "nAtom", nAtom_);
      #ifdef SIMP_BOND
      loadParameter<int>(ar, "nBond", nBond_);
      #endif
      #ifdef SIMP_ANGLE
      loadParameter<int>(ar, "nAngle", nAngle_);
      #endif
      #ifdef SIMP_DIHEDRAL
      loadParameter<int>(ar, "nDihedral", nDihedral_);
      #endif
      allocate();

      loadDArray<int>(ar, "atomTypeIds", atomTypeIds_, nAtom_);

      #ifdef SIMP_BOND
      if (nBond_ > 0) {
         loadDArray<SpeciesBond>(ar, "speciesBonds", speciesBonds_, 
                                 nBond_);
         // Make atomBondIdArrays
         int bondId, atomId1, atomId2;
         for (bondId = 0; bondId < nBond_; ++bondId) {
            atomId1 = speciesBond(bondId).atomId(0);
            atomId2 = speciesBond(bondId).atomId(1);
            atomBondIdArrays_[atomId1].append(bondId);
            atomBondIdArrays_[atomId2].append(bondId);
         }

      }
      #endif

      #ifdef SIMP_ANGLE
      if (nAngle_ > 0) {
         loadDArray<SpeciesAngle>(ar, "speciesAngles", speciesAngles_,
                                  nAngle_);
         // Make atomAngleIdArrays
         int angleId, atomId1, atomId2, atomId3;
         for (angleId = 0; angleId < nAngle_; ++angleId) {
            atomId1 = speciesAngle(angleId).atomId(0);
            atomId2 = speciesAngle(angleId).atomId(1);
            atomId3 = speciesAngle(angleId).atomId(2);
            atomAngleIdArrays_[atomId1].append(angleId);
            atomAngleIdArrays_[atomId2].append(angleId);
            atomAngleIdArrays_[atomId3].append(angleId);
         }
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      // Make atomDihedralIdArrays
      if (nDihedral_ > 0) {
         loadDArray<SpeciesDihedral>(ar, "speciesDihedrals", speciesDihedrals_,
                                     nDihedral_);
         int dihedralId;
         int tAtomId1, tAtomId2, tAtomId3, tAtomId4;
         for (dihedralId = 0; dihedralId < nDihedral_; ++dihedralId) {
            tAtomId1 = speciesDihedral(dihedralId).atomId(0);
            tAtomId2 = speciesDihedral(dihedralId).atomId(1);
            tAtomId3 = speciesDihedral(dihedralId).atomId(2);
            tAtomId4 = speciesDihedral(dihedralId).atomId(3);
            atomDihedralIdArrays_[tAtomId1].append(dihedralId);
            atomDihedralIdArrays_[tAtomId2].append(dihedralId);
            atomDihedralIdArrays_[tAtomId3].append(dihedralId);
            atomDihedralIdArrays_[tAtomId4].append(dihedralId);
         }
      }
      #endif
   }

   /*
   * Save internal state to an archive.
   */
   void Species::save(Serializable::OArchive &ar)
   {
      ar << id_;
      ar << moleculeCapacity_;
      ar << nAtom_;
      #ifdef SIMP_BOND
      ar << nBond_;
      #endif
      #ifdef SIMP_ANGLE
      ar << nAngle_;
      #endif
      #ifdef SIMP_DIHEDRAL
      ar << nDihedral_;
      #endif
      ar << atomTypeIds_;
      #ifdef SIMP_BOND
      ar << speciesBonds_;
      #endif
      #ifdef SIMP_ANGLE
      ar << speciesAngles_;
      #endif
      #ifdef SIMP_DIHEDRAL
      ar << speciesDihedrals_;
      #endif
   }

   // Setters

   /*
   * Set value of integer id for this species.
   */
   void Species::setId(int id)
   { 
      // Preconditions
      if (id < 0)   UTIL_THROW("Negative species id");
      if (id_ >= 0) UTIL_THROW("Species id was set previously");

      id_ = id; 
   }

   /*
   * Set a pointer to an associated SpeciesMutator object.
   */
   void Species::setMutatorPtr(McMd::SpeciesMutator* mutatorPtr)
   {  mutatorPtr_ = mutatorPtr; }

   /*
   * Allocate memory for arrays that describe chemical structure.
   */
   void Species::allocate() 
   {
      assert(nAtom_ >  0);
      atomTypeIds_.allocate(nAtom_);

      #ifdef SIMP_BOND
      atomBondIdArrays_.allocate(nAtom_);
      assert(nBond_ >= 0);
      if (nBond_ > 0) {
         speciesBonds_.allocate(nBond_);
      } 
      #endif
      #ifdef SIMP_ANGLE
      atomAngleIdArrays_.allocate(nAtom_);
      assert(nAngle_ >= 0);
      if (nAngle_ > 0) {
         speciesAngles_.allocate(nAngle_);
      } 
      #endif
      #ifdef SIMP_DIHEDRAL
      atomDihedralIdArrays_.allocate(nAtom_);
      assert(nDihedral_ >= 0);
      if (nDihedral_ > 0) {
         speciesDihedrals_.allocate(nDihedral_);
      } 
      #endif

      // Initialize atom type Ids to null/invalid values
      for (int i = 0; i < nAtom_; ++i) {
         atomTypeIds_[i] = -1;
      }
   }

   /*
   * Set the atom type for one atom
   */
   void Species::setAtomType(int atomId, int atomType)
   {  
      atomTypeIds_[atomId] = atomType;
   }

   #ifdef SIMP_BOND
   /*
   * Add a bond to the species chemical structure.
   */
   void Species::makeBond(int bondId, int atomId1, int atomId2, int bondType)
   {
      // Preconditions
      assert(bondId  >= 0);
      assert(atomId1 >= 0);
      assert(atomId2 >= 0);
      assert(bondId  < nBond_);
      assert(atomId1 < nAtom_);
      assert(atomId2 < nAtom_);

      // Create a SpeciesGroup<2> object for this bond
      speciesBonds_[bondId].setAtomId(0, atomId1);
      speciesBonds_[bondId].setAtomId(1, atomId2); 
      speciesBonds_[bondId].setTypeId(bondType); 

      // Add bond Id to AtomBondIdArray objects for both atoms
      atomBondIdArrays_[atomId1].append(bondId);
      atomBondIdArrays_[atomId2].append(bondId);

   }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Add an angle to the species chemical structure.
   */
   void Species::makeAngle(
       int angleId, int atomId1, int atomId2, int atomId3, int angleType)
   {
      // Preconditions.
      assert(angleId >= 0);
      assert(atomId1 >= 0);
      assert(atomId2 >= 0);
      assert(atomId3 >= 0);
      assert(angleId < nAngle_);
      assert(atomId1 < nAtom_);
      assert(atomId2 < nAtom_);
      assert(atomId3 < nAtom_);

      // Create a SpeciesGroup<3> object for this angle.
      speciesAngles_[angleId].setAtomId(0, atomId1);
      speciesAngles_[angleId].setAtomId(1, atomId2);
      speciesAngles_[angleId].setAtomId(2, atomId3);
      speciesAngles_[angleId].setTypeId(angleType);

      // Add angle Id to AtomAngleIdArray objects for all atoms.
      atomAngleIdArrays_[atomId1].append(angleId);
      atomAngleIdArrays_[atomId2].append(angleId);
      atomAngleIdArrays_[atomId3].append(angleId);
   }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Add a dihedral to the species chemical structure.
   */
   void Species::makeDihedral(int dihedralId, int atomId1, int atomId2,
                             int atomId3, int atomId4, int dihedralType)
   {
      // Preconditions.
      assert(dihedralId >= 0);
      assert(atomId1 >= 0);
      assert(atomId2 >= 0);
      assert(atomId3 >= 0);
      assert(atomId4 >= 0);
      assert(dihedralId < nDihedral_);
      assert(atomId1 < nAtom_);
      assert(atomId2 < nAtom_);
      assert(atomId3 < nAtom_);
      assert(atomId4 < nAtom_);

      // Create a SpeciesGroup<4> object for this angle.
      speciesDihedrals_[dihedralId].setAtomId(0, atomId1);
      speciesDihedrals_[dihedralId].setAtomId(1, atomId2);
      speciesDihedrals_[dihedralId].setAtomId(2, atomId3);
      speciesDihedrals_[dihedralId].setAtomId(3, atomId4);
      speciesDihedrals_[dihedralId].setTypeId(dihedralType);

      // Add angle Id to AtomAngleIdArray objects for all atoms.
      atomDihedralIdArrays_[atomId1].append(dihedralId);
      atomDihedralIdArrays_[atomId2].append(dihedralId);
      atomDihedralIdArrays_[atomId3].append(dihedralId);
      atomDihedralIdArrays_[atomId4].append(dihedralId);
   }
   #endif

   /*
   * Check validity of Species. Return true if valid, or throw an Exception.
   */
   bool Species::isValid() const
   {

      if (atomTypeIds_.isAllocated()) {
 
         // Check atomTypeIds array
         if (!isMutable()) {
            for (int i = 0; i < nAtom_; ++i) {
               if (atomTypeIds_[i] < 0) {
                  UTIL_THROW("Negative AtomTypeIds_ in non-mutable Species");
               }
            }
         }

         #ifdef SIMP_BOND
         {
            // Loop over all bonds (if any) in speciesBonds_ array
            int  atomId, bondId, atomId0, atomId1, j;
            bool hasBond;
            for (bondId = 0; bondId < nBond_; ++bondId) {
   
               atomId0 = speciesBond(bondId).atomId(0);
               atomId1 = speciesBond(bondId).atomId(1);
      
               // Check that atomIds are in valid range and unequal
               if (atomId0 < 0 || atomId0 >= nAtom_) {
                  UTIL_THROW("Invalid atomId0 in speciesBonds_");
               }
               if (atomId1 < 0 || atomId1 >= nAtom_) {
                  UTIL_THROW("Invalid atomId1 in speciesBonds_");
               }
               if (atomId0 == atomId1) {
                  UTIL_THROW("Equal atom ids in a SpeciesBond");
               }
   
               // Check that bondId is in atombondIdArrays_ for atomId0
               hasBond = false;
               for (j = 0; j < atomBondIdArrays_[atomId0].size(); ++j) {
                  if (atomBondIdArrays_[atomId0][j] == bondId) {
                     hasBond = true;
                     break;
                  } 
               }
               if (!hasBond) {
                  UTIL_THROW("BondId missing from atomBondIdArrays_");
               }
   
               // Check that bondId is in atombondIdArrays_ for atomId1
               hasBond = false;
               for (j = 0; j < atomBondIdArrays_[atomId1].size(); ++j) {
                  if (atomBondIdArrays_[atomId1][j] == bondId) {
                     hasBond = true;
                     break;
                  } 
               }
               if (!hasBond) {
                  UTIL_THROW("BondId missing from atomBondIdArrays_");
               }
   
            }
   
            // Loop over atomBondIdArrays for all atoms
            int bondCounter = 0;
            for (atomId = 0; atomId < nAtom_; ++atomId) {
               for (j = 0; j < atomBondIdArrays_[atomId].size(); ++j) {
                  bondId = atomBondIdArrays_[atomId][j];
   
                  if (bondId < 0 || bondId >= nBond_) 
                     UTIL_THROW("Bond index out of range");
   
                  atomId0 = speciesBond(bondId).atomId(0);
                  atomId1 = speciesBond(bondId).atomId(1);
                  if (atomId != atomId0 && atomId != atomId1) 
                     UTIL_THROW("Inconsistency in atomBondId and bond");
                  ++bondCounter;
                  
               }
            }
            if (bondCounter != 2*nBond_) {
               UTIL_THROW("Inconsistency in total number of bonds");
            }
         }
         #endif

         #ifdef SIMP_ANGLE
         {
            // Loop over all angles (if any) in speciesAngles_ array
            int  angleId, id, id2, atomId, atomId2(-1), j;
            bool hasAngles;
            for (angleId = 0; angleId < nAngle_; ++angleId) {
   
               for (id = 0; id < 3; ++id) {
                  atomId = speciesAngle(angleId).atomId(id);
   
                  // Check that atomIds are in valid range and unequal
                  if (atomId < 0 || atomId >= nAtom_) {
                     UTIL_THROW("Invalid atomId in speciesAngles_");
                  }
                  for (id2 = id + 1; id2 < 3; ++id2) {
                     atomId2 = speciesAngle(angleId).atomId(id2);
                     if (atomId2 == atomId) {
                         UTIL_THROW("Equal atom ids in a SpeciesAngle");
                     }
                  }
   
                  // Check that angleId is in atomAngleIdArrays_ for atomId
                  hasAngles = false;
                  for (j = 0; j < atomAngleIdArrays_[atomId].size(); ++j) {
                     if (atomAngleIdArrays_[atomId][j] == angleId) {
                        hasAngles = true;
                        break;
                     } 
                  }
                  if (!hasAngles) {
                     UTIL_THROW("AngleId missing from atomAngleIdArrays_");
                  }
               }
   
            }
   
            // Loop over atomAngleIdArrays for all atoms
            int angleCounter = 0;
            bool hasAtom;
            for (atomId = 0; atomId < nAtom_; ++atomId) {
               for (j = 0; j < atomAngleIdArrays_[atomId].size(); ++j) {
                  angleId = atomAngleIdArrays_[atomId][j];
                  if (angleId < 0 || angleId >= nAngle_) 
                     UTIL_THROW("Angle index out of range");
   
                  hasAtom = false;
                  for (id = 0; id < 3; ++id) {
                     atomId2 = speciesAngle(angleId).atomId(id);
                     if (atomId2 == atomId) {
                        hasAtom = true;
                        break;
                     }
                  }
                  if (!hasAtom) {
                     UTIL_THROW("Inconsistency in atomAngleId and angle");
                  }
                  ++angleCounter;
               }
            }
            if (angleCounter != 3*nAngle_) {
               UTIL_THROW("Inconsistency in total number of angles");
            }
         }
         #endif

         #ifdef SIMP_DIHEDRAL
         {
            // Loop over all dihedrals (if any) in speciesDihedrals_ array
            int  dihedralId, tId, tId2, tAtomId, tAtomId2, j;
            bool hasDihedral;
   
            for (dihedralId = 0; dihedralId < nDihedral_; ++dihedralId) {
   
               for (tId = 0; tId < 4; ++tId) {
                  tAtomId = speciesDihedral(dihedralId).atomId(tId);
   
                  // Check that atomIds are in valid range and unequal
                  if (tAtomId < 0 || tAtomId >= nAtom_) {
                     UTIL_THROW("Invalid atomId in speciesDihedral_");
                  }
                  for (tId2 = tId + 1; tId2 < 4; ++tId2) {
                     tAtomId2 = speciesDihedral(dihedralId).atomId(tId2);
                     if (tAtomId2 == tAtomId) {
                         UTIL_THROW("Equal atom ids in a SpeciesDihedral");
                     }
                  }
   
                  // Check that dihedralId is in atomDihedralIdArrays_ for tAtomId
                  hasDihedral = false;
                  for (j = 0; j < atomDihedralIdArrays_[tAtomId].size(); ++j) {
                     if (atomDihedralIdArrays_[tAtomId][j] == dihedralId) {
                        hasDihedral = true;
                        break;
                     } 
                  }
                  if (!hasDihedral) {
                     UTIL_THROW("DihedralId missing from atomDihedralIdArrays_");
                  }
               }
   
            }
   
            // Loop over atomDihedralIdArrays for all atoms.
            int dihedralCounter = 0;
            bool tHasAtom;
            for (tAtomId = 0; tAtomId < nAtom_; ++tAtomId) {
               for (j = 0; j < atomDihedralIdArrays_[tAtomId].size(); ++j) {
                  dihedralId = atomDihedralIdArrays_[tAtomId][j];
                  if (dihedralId < 0 || dihedralId >= nDihedral_) 
                     UTIL_THROW("Dihedral index out of range");
   
                  tHasAtom = false;
                  for (tId = 0; tId < 4; ++tId) {
                     tAtomId2 = speciesDihedral(dihedralId).atomId(tId);
                     if (tAtomId2 == tAtomId) {
                        tHasAtom = true;
                        break;
                     }
                  }
                  if (!tHasAtom) {
                     UTIL_THROW("Inconsistency in atomDihedralId and dihedral");
                  }
                  ++dihedralCounter;
               }
            }
            if (dihedralCounter != 4*nDihedral_) {
               UTIL_THROW("Inconsistency in total number of dihedrals");
            }
         }
         #endif

      } else { 

         if (moleculeCapacity_ != 0) {
            UTIL_THROW("Species not allocated but moleculeCapacity != 0");
         }

      }

      return true;
   }

} 
