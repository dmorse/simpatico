/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesStorage.h"
#include <mdPp/chemistry/Atom.h>

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   SpeciesStorage::SpeciesStorage()
    : atomPtrs_(),
      molecules_()
   {}

   /*
   * Destructor.
   */
   SpeciesStorage::~SpeciesStorage()
   {}

   void SpeciesStorage::initialize()
   {
      // Preconditions
      if (capacity() <= 0) {
         UTIL_THROW("SpeciesStorage capacity() not yet set");
      }
      if (molecules_.capacity() != 0) {
         UTIL_THROW("SpeciesStorage::molecules_ already allocated");
      }
      if (atomPtrs_.capacity() != 0) {
         UTIL_THROW("SpeciesStorage::atomPtrs_ already allocated");
      }
      if (capacity() <= 0) {
         UTIL_THROW("SpeciesStorage::capacity() <= 0: Value not set");
      }
      if (nAtom() <= 0) {
         UTIL_THROW("SpeciesStorage::nAtom() <= 0: Value not set");
      }

      // Allocate and initialize atomPtrs_ array
      atomPtrs_.allocate(capacity()*nAtom());
      for (int i=0; i < atomPtrs_.capacity(); ++i) {
         atomPtrs_[i] = 0;
      }

      // Allocate and initialize molecules_ array
      molecules_.allocate(capacity());
      molecules_.resize(capacity()); // temporarily resize
      Atom** atomPtr = &atomPtrs_[0];
      for (int i = 0; i < capacity(); ++i) {
         molecules_[i].atoms_ = atomPtr;
         molecules_[i].id_ = i;
         molecules_[i].nAtom_ = 0;
         molecules_[i].speciesPtr_ = this;
         atomPtr += nAtom();
      }
      molecules_.resize(0); // reset size to zero

   }

   /*
   * Reset to empty state.
   */
   void SpeciesStorage::clear()
   {
      for (int i = 0; i < atomPtrs_.capacity(); ++i) {
         atomPtrs_[i] = 0;
      }
      molecules_.resize(capacity()); // temporarily resize
      for (int i = 0; i < molecules_.capacity(); ++i) {
         molecules_[i].nAtom_ = 0;
      }
      molecules_.resize(0);
   }

   /*
   * Add an atom to this SpeciesStorage.
   */
   void SpeciesStorage::addAtom(Atom& atom)
   {
      int sId = atom.speciesId;
      if (sId != id_) {
         UTIL_THROW("Inconsistent speciesId");
      }
      int mId = atom.moleculeId;
      int aId = atom.atomId;
      if (mId < 0) {
         UTIL_THROW("atom.moleculeId < 0");
      }
      if (mId >= capacity()) {
         UTIL_THROW("atom.moleculeId >= capacity()");
      }
      if (aId < 0) {
         UTIL_THROW("atom.atomId < 0");
      }
      if (aId >= nAtom_) {
         UTIL_THROW("atom.atomId >= nAtom_");
      }
      int gid = mId*nAtom_ + aId; // "global" atom id, within species
      if (atomPtrs_[gid] != 0) {
         UTIL_THROW("Atom already present");
      }
      atomPtrs_[gid] = &atom;
      if (mId >= molecules_.size()) {
         molecules_.resize(mId+1);
      }
      ++molecules_[mId].nAtom_;
   }

   /*
   * Initialized an iterator over molecules in species.
   */
   void SpeciesStorage::begin(Iterator& iterator)
   {  molecules_.begin(iterator); }

   /*
   * Return true if valid, throw Exception otherwise.
   */
   bool SpeciesStorage::isValid() const
   {
      // Check allocation and initialization
      UTIL_CHECK(capacity() > 0);
      if (molecules_.capacity() != capacity()) {
         UTIL_THROW("molecules_ is not allocated");
      }
      if (atomPtrs_.capacity() != capacity()*nAtom()) {
         UTIL_THROW("atomPtrs_ is not allocated");
      }

      // Check molecules and atoms
      if (molecules_.size() > 0) {
         const Molecule* mPtr;
         const Atom* aPtr;
         int im, ia;
         for (im = 0; im < molecules_.size(); ++im) {
            mPtr = &(molecules_[im]);
            if (mPtr->atoms_ != &atomPtrs_[0] + im*nAtom_) {
               UTIL_THROW("Incorrect assignment of atoms_ in molecule");
            }
            if (mPtr->id() != im) {
               UTIL_THROW("Incorrect molecule id");
            }
            if (mPtr->speciesPtr_ != this) {
               UTIL_THROW("Incorrect species pointer in molecule");
            }
            if (mPtr->nAtom_ == nAtom_) {
               for (ia = 0; ia < nAtom_; ++ia) {
                  aPtr = mPtr->atoms_[ia];
                  if (aPtr == 0) {
                     UTIL_THROW("Null atom ptr in molecule");
                  }
                  if (aPtr->atomId != ia) {
                     UTIL_THROW("Inconsistent atom index");
                  }
                  if (aPtr->moleculeId != im) {
                     UTIL_THROW("Inconsistent molecule index");
                  }
                  if (aPtr->speciesId != id_) {
                     UTIL_THROW("Inconsistent species index");
                  }
               } 
            } else 
            if (mPtr->nAtom_ == 0) {
               for (ia = 0; ia < nAtom_; ++ia) {
                  if (mPtr->atoms_[ia] != 0) {
                     UTIL_THROW("Nonnull atom ptr in empty molecule");
                  }
               }
            } else {
               UTIL_THROW("0 != molecule nAtom != species nAtom");
            }
         }
      }
      return true;
   }

}
