/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Configuration.h"

namespace MdPp 
{

   /*
   * Constructor.
   */
   Configuration::Configuration()
    : 
      nSpecies_(0),
      atomCapacity_(0)
      #ifdef SIMP_BOND
      , bondCapacity_(0)
      #endif
      #ifdef SIMP_ANGLE
      , angleCapacity_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralCapacity_(0)
      #endif
      #ifdef SIMP_IMPROPERT
      , improperCapacity_(0)
      #endif
      , hasAtomContexts_(false)
   {  setClassName("Configuration"); }

   /*
   * Destructor.
   */
   Configuration::~Configuration()
   {}

   /*
   * Open, read and close a parameter file.
   */
   void Configuration::readParam(const char* filename)
   {
      std::ifstream in;
      in.open(filename);
      readParam(in);
      in.close();
   }

   /*
   * Read parameters from file.
   */
   void Configuration::readParameters(std::istream& in)
   {

      // Optionally read atom capacity
      atomCapacity_ = 0; // default value
      readOptional<int>(in, "atomCapacity", atomCapacity_); 
      if (atomCapacity_ > 0) {
         atoms_.allocate(atomCapacity_);
      }

      #ifdef SIMP_BOND
      bondCapacity_ = 0; // default value
      readOptional<int>(in, "bondCapacity", bondCapacity_); 
      // If bondCapacity is absent, it is set to zero by default
      if (bondCapacity_ > 0) {
         bonds_.allocate(bondCapacity_);
      }
      #endif

      #ifdef SIMP_ANGLE
      angleCapacity_ = 0; // default value
      readOptional<int>(in, "angleCapacity", angleCapacity_); 
      if (angleCapacity_ > 0) {
         angles_.allocate(angleCapacity_);
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      dihedralCapacity_ = 0; // default value
      readOptional<int>(in, "dihedralCapacity", dihedralCapacity_); 
      if (dihedralCapacity_ > 0) {
         dihedrals_.allocate(dihedralCapacity_);
      }
      #endif

      #ifdef SIMP_IMPROPER
      improperCapacity_ = 0; // default value
      readOptional<int>(in, "improperCapacity", improperCapacity_); 
      if (improperCapacity_ > 0) {
         impropers_.allocate(improperCapacity_);
      }
      #endif

   }

   void Configuration::setNSpecies(int nSpecies)
   {
      UTIL_CHECK(nSpecies > 0);
      nSpecies_ = nSpecies;
      species_.allocate(nSpecies_);
      for (int i = 0; i < nSpecies_; ++i) {
         species_[i].setId(i);
      }
   }

   /*
   * Set value of hasAtomContexts flag.
   */
   void Configuration::setHasAtomContexts(bool hasAtomContexts)
   {  hasAtomContexts_ = hasAtomContexts; }

   /*
   * Remove all atoms and groups - set to empty state.
   */
   void Configuration::clear()
   {
      // Clear species data
      if (nSpecies_ > 0) {
         UTIL_CHECK(species_.isAllocated());
         UTIL_CHECK(nSpecies_ = species_.capacity());
         for (int i = 0; i < nSpecies_; ++i) {
            species(i).clear();
         }
      } else {
         UTIL_CHECK(!species_.isAllocated());
      }

      // Clear atom data
      if (atoms_.capacity() > 0) {
         atoms_.clear();
      }
      hasAtomContexts_ = false;

      // Clear group data
      #ifdef SIMP_BOND
      if (bonds_.capacity() > 0) {
         bonds_.clear();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (angles_.capacity() > 0) {
         angles_.clear();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (dihedrals_.capacity() > 0) {
         dihedrals_.clear();
      }
      #endif
      #ifdef SIMP_IMPROPER
      if (impropers_.capacity() > 0) {
         impropers_.clear();
      }
      #endif
   }

   /*
   * Set atom context info, assuming ordered atom ids.
   */
   void Configuration::setAtomContexts()
   {
      UTIL_CHECK(nSpecies() > 0);

      // Check total number of atoms
      int nFill = 0;
      for (int i = 0; i < nSpecies(); ++i) {
         nFill += species(i).nAtom()*species(i).capacity();
      }
      if (atoms().size() != nFill) {
         UTIL_THROW("Atom storage capacity != expected total");
      }

      // Set speciesId, moleculeId and atomId for each Atom
      Atom* atomPtr = 0;
      SpeciesStorage* speciesPtr = 0;
      int speciesId, moleculeId, atomId, nMol, nAtom;
      int id = 0;
      for (speciesId = 0; speciesId < nSpecies(); ++speciesId) {
         speciesPtr = &species(speciesId);
         nMol  = speciesPtr->capacity();
         nAtom = speciesPtr->nAtom();
         for (moleculeId = 0; moleculeId < nMol; ++moleculeId) {
            for (atomId = 0; atomId < nAtom; ++atomId) {
               atomPtr = atoms().ptr(id);
               if (!atomPtr) {
                  UTIL_THROW("Error: missing atom");
               }
               atomPtr->speciesId = speciesId;
               atomPtr->moleculeId = moleculeId;
               atomPtr->atomId = atomId;
               ++id;
            }
         }
      }

      hasAtomContexts_ = true;
   }

   /*
   * Add all atoms currently in the AtomStorage to the Species 
   * containers and associated molecules.
   */
   void Configuration::addAtomsToSpecies()
   {
      UTIL_CHECK(nSpecies() > 0);
      for (int i = 0; i < nSpecies(); ++i) {
         species(i).clear();
      }
      int speciesId;
      AtomStorage::Iterator iter;
      atoms().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         speciesId = iter->speciesId;
         if (speciesId < 0) {
            UTIL_THROW("Negative speciesId");
         }
         if (speciesId >= nSpecies()) {
            UTIL_THROW("speciesId >= nSpecies");
         }
         species(speciesId).addAtom(*iter);
      }
   }


   void Configuration::makeGroups()
   {
      #ifdef SIMP_BOND
      makeBonds();
      #endif
      #ifdef SIMP_ANGLE
      makeAngles();
      #endif
      #ifdef SIMP_DIHEDRAL
      makeDihedrals();
      #endif
   }

   #ifdef SIMP_BOND
   /*
   * Create all bonds from species templates.
   */
   void Configuration::makeBonds()
   {
      UTIL_CHECK(nSpecies());

      // Count the required number of bonds
      SpeciesStorage* speciesPtr = 0;
      int nMol, nGroup;
      int nGroupTot = 0;
      for (int speciesId = 0; speciesId < nSpecies(); ++speciesId) {
         speciesPtr = &species(speciesId);
         nMol = speciesPtr->capacity();
         nGroup = speciesPtr->nBond();
         nGroupTot += nMol*nGroup;
      }
      if (nGroupTot == 0) return;

      // Allocate if necessary
      GroupStorage<2>& storage = bonds();
      if (storage.capacity() == 0) {
         storage.allocate(nGroupTot);
      }
      UTIL_CHECK(storage.capacity() >= nGroupTot);

      // Make groups
      int firstAtomId = 0;
      int groupId = 0;
      int nAtom;
      for (int speciesId = 0; speciesId < nSpecies(); ++speciesId) {
         speciesPtr = &species(speciesId);
         nMol = speciesPtr->capacity();
         nGroup = speciesPtr->nBond();
         nAtom = speciesPtr->nAtom();
         makeSpeciesGroups<2>(storage, speciesPtr->speciesBonds(),
                              nMol, nAtom, nGroup, firstAtomId, groupId);
      }
   }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Create all angles from species templates.
   */
   void Configuration::makeAngles()
   {
      UTIL_CHECK(nSpecies() > 0);

      // Count the required number of angles
      SpeciesStorage* speciesPtr = 0;
      int nMol, nGroup, nAtom;
      int nGroupTot = 0;
      for (int speciesId = 0; speciesId < nSpecies(); ++speciesId) {
         speciesPtr = &species(speciesId);
         int nMol = speciesPtr->capacity();
         int nGroup = speciesPtr->nAngle();
         nGroupTot += nMol*nGroup;
      }
      if (nGroupTot == 0) return;

      // Allocate if necessary
      GroupStorage<3>& storage = angles();
      if (storage.capacity() == 0) {
         storage.allocate(nGroupTot);
      }
      UTIL_CHECK(storage.capacity() >= nGroupTot);

      // Make groups
      int firstAtomId = 0;
      int groupId = 0;
      for (int speciesId = 0; speciesId < nSpecies(); ++speciesId) {
         speciesPtr = &species(speciesId);
         nMol = speciesPtr->capacity();
         nGroup = speciesPtr->nAngle();
         nAtom = speciesPtr->nAtom();
         makeSpeciesGroups<3>(storage, speciesPtr->speciesAngles(),
                              nMol, nAtom, nGroup, firstAtomId, groupId);
      }
   }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Create all dihedrals from species templates.
   */
   void Configuration::makeDihedrals()
   {
      UTIL_CHECK(nSpecies());

      // Count the required number of dihedrals
      SpeciesStorage* speciesPtr = 0;
      int nMol, nGroup, nAtom;
      int nGroupTot = 0;
      for (int speciesId = 0; speciesId < nSpecies(); ++speciesId) {
         speciesPtr = &species(speciesId);
         int nMol = speciesPtr->capacity();
         int nGroup = speciesPtr->nDihedral();
         nGroupTot += nMol*nGroup;
      }
      if (nGroupTot == 0) return;

      // Allocate if necessary
      GroupStorage<4>& storage = dihedrals();
      if (storage.capacity() == 0) {
         storage.allocate(nGroupTot);
      }
      UTIL_CHECK(storage.capacity() >= nGroupTot);

      // Make groups
      int firstAtomId = 0;
      int groupId = 0;
      for (int speciesId = 0; speciesId < nSpecies(); ++speciesId) {
         speciesPtr = &species(speciesId);
         nMol = speciesPtr->capacity();
         nGroup = speciesPtr->nDihedral();
         nAtom = speciesPtr->nAtom();
         makeSpeciesGroups<4>(storage, speciesPtr->speciesDihedrals(),
                              nMol, nAtom, nGroup, firstAtomId, groupId);
      }
   }
   #endif

   /**
   * Make all Group<N> objects of one type for one species.
   */
   template <int N>
   void Configuration::makeSpeciesGroups(
       GroupStorage<N>& storage,
       const DArray< SpeciesGroup<N> >& speciesGroups,
       int nMol, int nAtom, int nGroup, 
       int& firstAtomId, int& groupId)
   {
      Group<N>* groupPtr = 0;
      const SpeciesGroup<N>* speciesGroupPtr = 0;
      int molId, atomId, typeId, i, j;

      for (molId = 0; molId < nMol; ++molId) {
         for (i = 0; i < nGroup; ++i) {
            groupPtr = storage.newPtr();
            groupPtr->id = groupId;
            speciesGroupPtr = &speciesGroups[i]; //
            typeId = speciesGroupPtr->typeId();
            groupPtr->typeId = typeId;
            for (j = 0; j < N; ++j) {
               atomId = firstAtomId + speciesGroupPtr->atomId(j);
               groupPtr->atomIds[j] = atomId;
            }
            ++groupId;
         }
         firstAtomId += nAtom;
      }
   }

}
