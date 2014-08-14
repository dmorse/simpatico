#ifndef SPAN_SPECIES_CPP
#define SPAN_SPECIES_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"
#include "Atom.h"

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   Species::Species()
    : atomPtrs_(),
      molecules_(),
      id_(-1),
      nAtom_(0),
      capacity_(0)
   {}

   void Species::setId(int id)
   { id_ = id; }

   void Species::initialize(int nAtom, int capacity)
   {
      if (capacity_ != 0) {
         UTIL_THROW("Species already initialized: capacity_ != 0");
      }
      if (nAtom_ != 0) {
         UTIL_THROW("Species already initialized: nAtom_ != 0");
      }
      nAtom_ = nAtom;
      capacity_ = capacity;
      initialize();
   }

   void Species::initialize()
   {
      // Preconditions
      if (molecules_.capacity() != 0) {
         UTIL_THROW("Species::molecules_ already allocated");
      }
      if (atomPtrs_.capacity() != 0) {
         UTIL_THROW("Species::atomPtrs_ already allocated");
      }
      if (capacity_ <= 0) {
         UTIL_THROW("Species::capacity_ <= 0: Value not set");
      }
      if (nAtom_ <= 0) {
         UTIL_THROW("Species::nAtom_ <= 0: Value not set");
      }

      // Allocate and initialize atomPtrs_ array
      atomPtrs_.allocate(capacity_*nAtom_);
      for (int i=0; i < atomPtrs_.capacity(); ++i) {
         atomPtrs_[i] = 0;
      }

      // Allocate and initialize molecules_ array
      molecules_.allocate(capacity_);
      molecules_.resize(capacity_); // temporarily resize
      Atom** atomPtr = &atomPtrs_[0];
      for (int i = 0; i < capacity_; ++i) {
         molecules_[i].atoms_ = atomPtr;
         molecules_[i].id_ = i;
         molecules_[i].nAtom_ = 0;
         molecules_[i].speciesPtr_ = this;
         atomPtr += nAtom_;
      }
      molecules_.resize(0); // reset size to zero

   }

   /*
   * Reset to empty state.
   */
   void Species::clear()
   {
      for (int i = 0; i < atomPtrs_.capacity(); ++i) {
         atomPtrs_[i] = 0;
      }
      for (int i = 0; i < molecules_.size(); ++i) {
         molecules_[i].nAtom_ = 0;
      }
      molecules_.resize(0);
   }

   /*
   * Add an atom to this Species.
   */
   void Species::addAtom(Atom& atom)
   {
      if (atom.speciesId != id_) {
         UTIL_THROW("Inconsistent speciesId");
      }
      int mId = atom.moleculeId;
      int aId = atom.atomId;
      if (mId < 0) {
         UTIL_THROW("atom.moleculeId < 0");
      }
      if (mId >= capacity_) {
         UTIL_THROW("atom.moleculeId >= capacity_");
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
   void Species::begin(MoleculeIterator& iterator)
   {  molecules_.begin(iterator); }

   /*
   * Return true if valid, throw Exception otherwise.
   */
   bool Species::isValid() const
   {
      const Molecule* mPtr;
      const Atom* aPtr;
      int ia, im;
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
      return true;
   }

   /*
   * istream extractor (>>) for a Species.
   */
   std::istream& operator >> (std::istream& in, Species &species)
   {
      in >> species.nAtom_;
      in >> species.capacity_;
      species.initialize();
      return in;
   }

   /*
   * ostream inserter (<<) for a Species.
   */
   std::ostream& operator << (std::ostream& out, const Species &species) 
   {
      out.width(10);
      out << species.nAtom_;
      out.width(10);
      out << species.capacity_;
      return out;
   }

}
#endif
