#ifndef DDMD_ATOM_STORAGE_CPP
#define DDMD_ATOM_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomStorage.h"
#include "AtomIterator.h"
#include "ConstAtomIterator.h"
#include "GhostIterator.h"
#include "ConstGhostIterator.h"
#include <ddMd/chemistry/Group.h>
#include <util/format/Int.h>
#include <util/mpi/MpiLoader.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   AtomStorage::AtomStorage()
    : atoms_(),
      atomSet_(),
      atomReservoir_(),
      ghosts_(),
      ghostSet_(),
      ghostReservoir_(),
      atomPtrs_(),
      newAtomPtr_(0),
      newGhostPtr_(0),
      atomCapacity_(0),
      ghostCapacity_(0),
      totalAtomCapacity_(0),
      maxNAtomLocal_(0),
      maxNGhostLocal_(0),
      #ifdef UTIL_MPI
      nAtomTotal_(),
      maxNAtom_(),
      maxNGhost_(),
      #endif
      locked_(false),
      isInitialized_(false),
      #if UTIL_ORTHOGONAL
      isCartesian_(true)
      #else
      isCartesian_(false)
      #endif
   {  setClassName("AtomStorage"); }
 
   /*
   * Destructor.
   */
   AtomStorage::~AtomStorage()
   {}

   /*
   * Set parameters and allocate memory.
   */
   void AtomStorage::initialize(int atomCapacity, int ghostCapacity, 
      int totalAtomCapacity)
   {
      atomCapacity_ = atomCapacity;
      ghostCapacity_ = ghostCapacity;
      totalAtomCapacity_ = totalAtomCapacity;
      allocate();
   }

   /*
   * Read parameters and allocate memory.
   */
   void AtomStorage::readParameters(std::istream& in)
   {
      read<int>(in, "atomCapacity", atomCapacity_);
      read<int>(in, "ghostCapacity", ghostCapacity_);
      read<int>(in, "totalAtomCapacity", totalAtomCapacity_);
      allocate();
   }

   /*
   * Load parameters and allocate memory (call on all processors)
   */
   void AtomStorage::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<int>(ar, "atomCapacity", atomCapacity_);
      loadParameter<int>(ar, "ghostCapacity", ghostCapacity_);
      loadParameter<int>(ar, "totalAtomCapacity", totalAtomCapacity_);
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(maxNAtomLocal_);
      loader.load(maxNGhostLocal_);
      allocate();
   }

   /*
   * Save parameters (call only on ioProcessor).
   */
   void AtomStorage::save(Serializable::OArchive& ar)
   {
      ar << atomCapacity_;
      ar << ghostCapacity_;
      ar << totalAtomCapacity_;
      ar << maxNAtomLocal_;
      ar << maxNGhostLocal_;
   }

   /*
   * Allocate and initialize all containers (private).
   */
   void AtomStorage::allocate()
   {
      // Precondition
      if (isInitialized_) {
         UTIL_THROW("AtomStorage can only be initialized once");
      }

      atoms_.allocate(atomCapacity_);
      atomReservoir_.allocate(atomCapacity_);
      atomSet_.allocate(atoms_);
      int i;
      for (i = atomCapacity_ - 1; i >=0; --i) {
          atomReservoir_.push(atoms_[i]);
      }

      ghosts_.allocate(ghostCapacity_);
      ghostReservoir_.allocate(ghostCapacity_);
      ghostSet_.allocate(ghosts_);
      for (i = ghostCapacity_ - 1; i >=0; --i) {
          ghostReservoir_.push(ghosts_[i]);
      }

      atomPtrs_.allocate(totalAtomCapacity_);
      for (i = 0; i < totalAtomCapacity_; ++i) {
         atomPtrs_[i] = 0;
      }

      snapshot_.allocate(atomCapacity_);

      isInitialized_ = true;
   }

   // Local atom mutators

   /*
   * Returns address for a new local Atom.
   */ 
   Atom* AtomStorage::newAtomPtr()
   {
      // Preconditions
      if (newAtomPtr_ != 0) 
         UTIL_THROW("Unregistered newAtomPtr_ still active");
      if (locked_ ) 
         UTIL_THROW("AtomStorage is locked");
      if (ghostSet_.size() > 0)
         UTIL_THROW("Ghosts are not cleared");
      if (atomReservoir_.size() == 0) 
         UTIL_THROW("Atom to pop from empty atomReservoir");

      newAtomPtr_ = &atomReservoir_.pop();
      newAtomPtr_->clear();
      return newAtomPtr_;
   }

   /*
   * Register new local Atom in internal data structures.
   */ 
   void AtomStorage::addNewAtom()
   {
      if (newAtomPtr_ == 0) 
         UTIL_THROW("No active newAtomPtr_");
      if (locked_) 
         UTIL_THROW("AtomStorage is locked");
      int atomId = newAtomPtr_->id();
      if (atomId < 0 || atomId >= totalAtomCapacity_) {
         std::cout << "atomId = " << atomId << std::endl;
         UTIL_THROW("atomId is out of range");
      }
      if (atomPtrs_[atomId] != 0) {
         std::cout << "atomId       = " << atomId << std::endl;
         std::cout << "New Position = " << newAtomPtr_->position() 
                   << std::endl;
         std::cout << "Old Position = " << atomPtrs_[atomId]->position() 
                   << std::endl;
         UTIL_THROW("Atom with specified id is already present");
      }

      // Add to atom set
      atomSet_.append(*newAtomPtr_);
      atomPtrs_[atomId] = newAtomPtr_;
      newAtomPtr_->setIsGhost(false);

      // De-activate new atom pointer
      newAtomPtr_ = 0;

      // Update statistics
      if (atomSet_.size() > maxNAtomLocal_) {
         maxNAtomLocal_ = atomSet_.size();
      }
   }

   Atom* AtomStorage::addAtom(int id)
   {
      Atom* ptr = newAtomPtr();
      ptr->setId(id);
      addNewAtom();
      return ptr;
   }

   /*
   * Remove a specific local Atom.
   */
   void AtomStorage::removeAtom(Atom* atomPtr)
   {
      // Precondition
      int atomId = atomPtr->id();
      if (atomId < 0 || atomId >= totalAtomCapacity_) {
         std::cout << "atomId = " << atomId << std::endl;
         UTIL_THROW("atomId is out of range");
      }
      atomReservoir_.push(*atomPtr);
      atomSet_.remove(*atomPtr);
      atomPtrs_[atomId] = 0;
      atomPtr->setId(-1);
      atomPtr->setIsGhost(false);
   }

   /*
   * Remove all local atoms.
   */
   void AtomStorage::clearAtoms()
   {
      // Precondition
      if (locked_) {
         UTIL_THROW("AtomStorage is locked");
      }

      Atom* atomPtr;
      int   atomId;
      while (atomSet_.size() > 0) {
         atomPtr = &atomSet_.pop();
         atomId = atomPtr->id();
         if (atomPtrs_[atomId] == atomPtr) { 
            atomPtrs_[atomId] = 0;
         } else {
            if (atomPtrs_[atomId] == 0) {
               UTIL_THROW("Error: Unexpected null in atomPtrs_");
            } else {
               if (atomPtrs_[atomId]->id() != atomId) {
                  UTIL_THROW("Error: Inconsistent id in atomPtrs_");
               }
            }
         }
         // atomPtr->setId(-1);
         // atomPtr->clear();
         atomReservoir_.push(*atomPtr);
      }

      // Sanity check.
      if (atomSet_.size() != 0) {
         UTIL_THROW("Nonzero atomSet size at end of clearAtoms");
      }
   }

   // Ghost atom mutators

   /*
   * Return address for a new ghost Atom.
   */ 
   Atom* AtomStorage::newGhostPtr()
   {
      // Preconditions
      if (newGhostPtr_ != 0) 
         UTIL_THROW("Unregistered newGhostPtr_ still active");
      if (locked_ ) 
         UTIL_THROW("AtomStorage is locked");
      if (ghostReservoir_.size() == 0) 
         UTIL_THROW("Atom to pop from empty atomReservoir");

      newGhostPtr_ = &ghostReservoir_.pop();
      newGhostPtr_->clear();
      newGhostPtr_->setIsGhost(true);
      return newGhostPtr_;
   }

   /*
   * Register new ghost Atom to internal data structures.
   */ 
   void AtomStorage::addNewGhost()
   {
      if (newGhostPtr_ == 0) {
         UTIL_THROW("No active newGhostPtr_");
      }
      if (locked_) {
         UTIL_THROW("AtomStorage is locked");
      }
      int atomId = newGhostPtr_->id();
      if (atomId < 0 || atomId >= totalAtomCapacity_) {
         UTIL_THROW("atomId is out of range");
      }

      // Add to ghost set
      ghostSet_.append(*newGhostPtr_);
      if (atomPtrs_[atomId] == 0) {
         atomPtrs_[atomId] = newGhostPtr_;
      }
      // Note another atom with same id may already be present.
      newGhostPtr_->setIsGhost(true);

      // De-activate new ghost pointer
      newGhostPtr_ = 0;

      // Update statistics
      if (ghostSet_.size() > maxNGhostLocal_) {
         maxNGhostLocal_ = ghostSet_.size();
      }
   }

   Atom* AtomStorage::addGhost(int id)
   {
      Atom* ptr = newGhostPtr();
      ptr->setId(id);
      addNewGhost();
      return ptr;
   }

   /*
   * Remove a specific ghost Atom.
   */
   void AtomStorage::removeGhost(Atom* atomPtr)
   {
      // Precondition
      if (locked_) 
         UTIL_THROW("AtomStorage is locked");

      ghostReservoir_.push(*atomPtr);
      ghostSet_.remove(*atomPtr);
      int atomId = atomPtr->id();
      if (atomPtrs_[atomId] == atomPtr) {
         atomPtrs_[atomId] = 0;
      }
      atomPtr->setId(-1);
   }

   /*
   * Remove all ghosts.
   */
   void AtomStorage::clearGhosts()
   {
      // Precondition
      if (locked_) 
         UTIL_THROW("AtomStorage is locked");

      Atom* atomPtr;
      int   atomId;
      while (ghostSet_.size() > 0) {
         atomPtr = &ghostSet_.pop();
         atomId = atomPtr->id();
         if (atomPtrs_[atomId] == atomPtr) { 
            atomPtrs_[atomId] = 0;
         } else {
            if (atomPtrs_[atomId] == 0) {
               UTIL_THROW("Error: Unexpected null in atomPtrs_");
            } else {
               if (atomPtrs_[atomId]->id() != atomId) {
                  UTIL_THROW("Error: Inconsistent id in atomPtrs_");
               }
            }
         }
         // atomPtr->setId(-1);
         // atomPtr->clear();
         ghostReservoir_.push(*atomPtr);
      }

      if (ghostSet_.size() != 0) {
         UTIL_THROW("Nonzero ghostSet size at end of clearGhosts");
      }
   }

   // Snapshot functions

   /*
   * Record current positions of all local atoms and lock storage.
   */
   void AtomStorage::makeSnapshot()
   {
      // Precondition
      if (!isCartesian()) {
         UTIL_THROW("Error: Coordinates not Cartesian in makeSnapshot");
      }
 
      AtomIterator iter;
      int i = 0;
      for (begin(iter); iter.notEnd(); ++iter) {
         snapshot_[i] = iter->position();
         ++i;
      }
      locked_ = true;
   }

   /*
   * Clear snapshot of local atom positions.
   */
   void AtomStorage::clearSnapshot()
   {
      locked_ = false; 
   }

   /*
   * Return max. sq. displacement of local atoms on this node since snapshot.
   */
   double AtomStorage::maxSqDisplacement()
   {
      if (!isCartesian()) {
         UTIL_THROW("Error: Coordinates not Cartesian in maxSqDisplacement");
      } 
      if (!locked_) {
         UTIL_THROW("Error: AtomStorage not locked in maxSqDisplacement");
      } 
      Vector dr;
      double norm;
      double max = 0.0;
      AtomIterator iter;
      int i = 0;
      for (begin(iter); iter.notEnd(); ++iter) {
         dr.subtract(iter->position(), snapshot_[i]);
         norm = dr.square();
         if (norm > max) {
            max = norm;
         }
         ++i;
      }
      return max;
   }

   // Accessors

   /*
   * Return pointer to Atom with specified id.
   */
   Atom* AtomStorage::find(int atomId) const
   {  return atomPtrs_[atomId]; }

   /*
   * Set iterator to beginning of the set of local atoms.
   */
   void AtomStorage::begin(AtomIterator& iterator)
   {  atomSet_.begin(iterator); }
 
   /*
   * Set const iterator to beginning of the set of local atoms.
   */
   void AtomStorage::begin(ConstAtomIterator& iterator) const
   {  atomSet_.begin(iterator); }
 
   /*
   * Set iterator to beginning of the set of ghost atoms.
   */
   void AtomStorage::begin(GhostIterator& iterator)
   {  ghostSet_.begin(iterator); }

   /*
   * Set iterator to beginning of the set of ghost atoms.
   */
   void AtomStorage::begin(ConstGhostIterator& iterator) const
   {  ghostSet_.begin(iterator); }

   /*
   * Transform all atomic coordinates from Cartesian to generalized.
   */
   void AtomStorage::transformCartToGen(const Boundary& boundary) 
   {
      if (!isCartesian()) {
         UTIL_THROW("Error: Coordinates not Cartesian on entry");
      }
      Vector r;
      if (nAtom()) {
         AtomIterator  atomIter;
         for (begin(atomIter); atomIter.notEnd(); ++atomIter) {
            r = atomIter->position();
            boundary.transformCartToGen(r, atomIter->position());
         }
      }
      if (nGhost()) {
         GhostIterator  ghostIter;
         for (begin(ghostIter); ghostIter.notEnd(); ++ghostIter) {
            r = ghostIter->position();
            boundary.transformCartToGen(r, ghostIter->position());
         }
      }
      isCartesian_ = false;
   }

   /*
   * Transform all atomic coordinates from generalized to Cartesian.
   */
   void AtomStorage::transformGenToCart(const Boundary& boundary) 
   {
      if (isCartesian()) {
         UTIL_THROW("Error: Coordinates are Cartesian on entry");
      }
      Vector r;
      if (nAtom()) {
         AtomIterator atomIter;
         for (begin(atomIter); atomIter.notEnd(); ++atomIter) {
            r = atomIter->position();
            boundary.transformGenToCart(r, atomIter->position());
         }
      }
      if (nGhost()) {
         GhostIterator ghostIter;
         for (begin(ghostIter); ghostIter.notEnd(); ++ghostIter) {
            r = ghostIter->position();
            boundary.transformGenToCart(r, ghostIter->position());
         }
      }
      isCartesian_ = true;
   }

   #ifdef UTIL_MPI
   /*
   * Compute, store and return total number of atoms on all processors.
   */
   void AtomStorage::computeNAtomTotal(MPI::Intracomm& communicator)
   {
      // If nAtomTotal is already set, do nothing and return.
      // if (nAtomTotal_.isSet()) return;

      int nAtomLocal = nAtom();
      int nAtomTotal = 0;
      communicator.Reduce(&nAtomLocal, &nAtomTotal, 1, 
                          MPI::INT, MPI::SUM, 0);
      if (communicator.Get_rank() !=0) {
         nAtomTotal = 0;
      }
      nAtomTotal_.set(nAtomTotal);
   }

   void AtomStorage::unsetNAtomTotal()
   {  nAtomTotal_.unset(); }
   #endif

   /*
   * Compute memory usage statistics (call on all processors).
   */
   #ifdef UTIL_MPI
   void AtomStorage::computeStatistics(MPI::Intracomm& communicator)
   #else
   void AtomStorage::computeStatistics()
   #endif
   { 
      #ifdef UTIL_MPI
      int maxNAtomGlobal;
      communicator.Allreduce(&maxNAtomLocal_, &maxNAtomGlobal, 1, 
                             MPI::INT, MPI::MAX);
      maxNAtom_.set(maxNAtomGlobal);
      maxNAtomLocal_ = maxNAtomGlobal;
      #else
      maxNAtom_.set(maxNAtomLocal_);
      #endif

      #ifdef UTIL_MPI
      int maxNGhostGlobal;
      communicator.Allreduce(&maxNGhostLocal_, &maxNGhostGlobal, 1, 
                             MPI::INT, MPI::MAX);
      maxNGhost_.set(maxNGhostGlobal);
      maxNGhostLocal_ = maxNGhostGlobal;
      #else
      maxNGhost_.set(maxNGhostLocal_);
      #endif
   }

   /*
   * Clear all statistics.
   */
   void AtomStorage::clearStatistics() 
   {
      maxNAtomLocal_ = 0;
      maxNAtom_.unset();
      maxNGhostLocal_ = 0;
      maxNGhost_.unset();
   }

   /*
   * Output statistics.
   */
   void AtomStorage::outputStatistics(std::ostream& out)
   {

      out << std::endl;
      out << "AtomStorage" << std::endl;
      out << "NAtom:  max, capacity     " 
                  << Int(maxNAtom_.value(), 10)
                  << Int(atomCapacity_, 10)
                  << std::endl;
      out << "NGhost: max, capacity     " 
                  << Int(maxNGhost_.value(), 10)
                  << Int(ghostCapacity_, 10)
                  << std::endl;
   }

   /*
   * Check validity of this AtomStorage.
   *
   * Returns true if all is ok, or throws an Exception.
   */
   bool AtomStorage::isValid() const
   {
      
      if (nAtom() + atomReservoir_.size() != atomCapacity_) 
         UTIL_THROW("nAtom + reservoir size != local capacity"); 

      if (nGhost() + ghostReservoir_.size() != ghostCapacity_) 
         UTIL_THROW("nGhost + reservoir size != ghost capacity"); 

      // Check consistency of indexing in atomPtrs_.
      // For atomPtr_[i] != 0, does atomPtr_[i]->id() = i ?
      Atom* ptr;
      int   i, j;
      j = 0;
      for (i = 0; i < totalAtomCapacity_ ; ++i) {
         ptr = atomPtrs_[i];
         if (ptr != 0) {
            ++j;
            if (ptr->id() != i) {
               std::cout << std::endl;
               std::cout << "Index i in atomPtrs_  " << i << std::endl;
               std::cout << "atomPtrs_[i]->id()    " << ptr->id() << std::endl;
               UTIL_THROW("ptr->id() != i");
            }
         }
      }

      // Iterate over, count and find local atoms on this processor.
      ConstAtomIterator localIter;
      j = 0;
      for (begin(localIter); localIter.notEnd(); ++localIter) {
         ++j;
         ptr = find(localIter->id());
         if (ptr == 0) {
            UTIL_THROW("Unable to find local atom returned by iterator"); 
         }
         if (ptr != localIter.get()) {
            UTIL_THROW("Inconsistent find(localIter->id()"); 
         }
      }
      if (j != nAtom()) {
         UTIL_THROW("Number counted by localIterator != nAtom()"); 
      }

      // Iterate over, count and find ghost atoms
      ConstGhostIterator ghostIter;
      j = 0;
      for (begin(ghostIter); ghostIter.notEnd(); ++ghostIter) {
         ++j;
         ptr = find(ghostIter->id());
         if (ptr == 0) {
            UTIL_THROW("find(ghostIter->id()) == 0"); 
         }
         // We do NOT test if ptr == ghostIter.get(), because it is possible
         // to have multiple atoms with the same id on one processor. One to 
         // one correspondence of ids and pointers is guaranteed only for 
         // local atoms.
      }
      if (j != nGhost()) {
         UTIL_THROW("Number counted by ghostIterator != nGhost()"); 
      }

      return true;
   }

   #ifdef UTIL_MPI
   /**
   * Return true if the internal state is valid, or throw an Exception.
   *  
   * \param communicator domain communicator for all domain processors.
   */
   bool AtomStorage::isValid(MPI::Intracomm& communicator) const
   {
      isValid();
      nAtomTotal_.isValid(communicator);
      return true;
   }
   #endif

}
#endif
