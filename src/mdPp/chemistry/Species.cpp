#ifndef MDPP_SPECIES_CPP
#define MDPP_SPECIES_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"
#include "Atom.h"

namespace MdPp
{

   using namespace Util;

   // Constructor.
   Species::Species()
    : atomPtrs_(),
      molecules_(),
      id_(-1),
      nAtom_(0),
      capacity_(0),
      size_(0)
   {}

   void Species::setId(int id)
   {  id_ = id; }

   void Species::initialize(int capacity, int nAtom)
   {
      capacity_ = capacity;
      nAtom_ = nAtom;
      initialize();
   }

   void Species::initialize()
   {
      // Preconditions
      if (molecules_.capacity() > 0) {
         UTIL_THROW("Species already initialized");
      }
      if (capacity_ <= 0) {
         UTIL_THROW("Capacity <= 0: Not initialized");
      }
      if (nAtom_ <= 0) {
         UTIL_THROW("nAtom_ <= 0: Not initialized");
      }

      // Allocate and initialize atomPtrs_ array
      atomPtrs_.allocate(capacity_*nAtom_);
      for (int i=0; i < atomPtrs_.capacity(); ++i) {
         atomPtrs_[i] = 0;
      }

      // Allocate and initialize molecules_ array
      molecules_.allocate(capacity_);
      Atom** atomPtr = &atomPtrs_[0];
      for (int i=0; i < capacity_; ++i) {
         molecules_[i].atoms_ = atomPtr;
         molecules_[i].id_ = i;
         molecules_[i].nAtom_ = 0;
         molecules_[i].speciesPtr_ = this;
         atomPtr += nAtom_;
      }
      size_ = 0;

   }

   void Species::clear()
   {
      for (int i=0; i < atomPtrs_.capacity(); ++i) {
         atomPtrs_[i] = 0;
      }
      for (int i=0; i < capacity_; ++i) {
         molecules_[i].nAtom_ = 0;
      }
      size_ = 0;
   }

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
      atomPtrs_[mId*nAtom_ + aId] = &atom;
      ++molecules_[mId].nAtom_;
      if (mId > size_) {
         size_ = mId;
      }
   }

   void Species::begin(MoleculeIterator& iterator)
   {  molecules_.begin(iterator); }

   bool Species::isValid() const
   {
      const Molecule* mPtr;
      const Atom* aPtr;
      int ia, im;
      for (im = 0; im < size_; ++im) {
         mPtr = &(molecules_[im]);
         if (mPtr->nAtom_ != nAtom_) {
            UTIL_THROW("molecule nAtom != species nAtom");
         }
         if (mPtr->atoms_ != &atomPtrs_[0] + im*nAtom_) {
            UTIL_THROW("Incorrect assignment of atoms_ in molecule");
         }
         if (mPtr->id() != im) {
            UTIL_THROW("Incorrect molecule id");
         }
         if (mPtr->speciesPtr_ != this) {
            UTIL_THROW("Incorrect species pointer in molecule");
         }
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
