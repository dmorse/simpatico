#ifndef DDMD_ATOM_STORAGE_CPP
#define DDMD_ATOM_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomStorage.h"
#include "AtomIterator.h"
#include "ConstAtomIterator.h"
#include "GhostIterator.h"
#include "ConstGhostIterator.h"
#include <ddMd/chemistry/Group.h>
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
      locked_(false),
      isInitialized_(false)
   {}
 
   /*
   * Destructor.
   */
   AtomStorage::~AtomStorage()
   {}

   /*
   * Read parameters and allocate memory.
   */
   void AtomStorage::readParam(std::istream& in)
   {
      readBegin(in, "AtomStorage");
      read<int>(in, "atomCapacity", atomCapacity_);
      read<int>(in, "ghostCapacity", ghostCapacity_);
      read<int>(in, "totalAtomCapacity", totalAtomCapacity_);
      allocate();
      readEnd(in);
   }

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

   /**
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

      atomSet_.append(*newAtomPtr_);
      atomPtrs_[atomId] = newAtomPtr_;
      newAtomPtr_->setIsGhost(false);
      newAtomPtr_ = 0;
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
      newGhostPtr_ = &ghostReservoir_.pop();
      newGhostPtr_->clear();
      return newGhostPtr_;
   }

   /*
   * Register new ghost Atom to internal data structures.
   */ 
   void AtomStorage::addNewGhost()
   {
      if (newGhostPtr_ == 0) 
         UTIL_THROW("No active newGhostPtr_");
      if (locked_) 
         UTIL_THROW("AtomStorage is locked");
      int atomId = newGhostPtr_->id();
      if (atomId < 0 || atomId >= totalAtomCapacity_) 
         UTIL_THROW("atomId is out of range");
      //if (atomPtrs_[atomId] != 0)
      //   UTIL_THROW("Atom with specified id is already present");

      ghostSet_.append(*newGhostPtr_);
      if (atomPtrs_[atomId] == 0) {
         atomPtrs_[atomId] = newGhostPtr_;
      }
      newGhostPtr_->setIsGhost(true);
      newGhostPtr_ = 0;
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
   * Remove a specific ghost Atom.
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
         }
         atomPtr->setId(-1);
         atomPtr->clear();
         ghostReservoir_.push(*atomPtr);
      }

      if (ghostSet_.size() != 0) {
         UTIL_THROW("Nonzero ghostSet size at end of clearGhosts");
      }
   }

   // Snapshot functions

   /**
   * Record current positions of all local atoms.
   */
   void AtomStorage::makeSnapshot()
   {
      AtomIterator iter;
      int i = 0;
      for (begin(iter); iter.notEnd(); ++iter) {
         snapshot_[i] = iter->position();
         ++i;
      }
      locked_ = true;
   }

   /**
   * Clear snapshot of local atom positions.
   */
   void AtomStorage::clearSnapshot()
   {  locked_ = false; }

   /**
   * Return maximum squared displacement on this processor since last snapshot.
   */
   double AtomStorage::maxSqDisplacement()
   {
      Vector dr;
      double norm;
      double max = 0;
      AtomIterator iter;
      int i = 0;
      for (begin(iter); iter.notEnd(); ++iter) {
         dr.subtract(iter->position(), snapshot_[i]);
         norm = dr.square();
         if (norm > max) max = norm;
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
               UTIL_THROW("ptr->id() != i"); 
            }
         }
      }
      // Feb. 2012: I don't remember why this test is commented out.
      // if (nGhost() + nAtom() != j) 
      //   UTIL_THROW("nGhost + nAtom != j"); 

      // Count local atoms on this processor.
      ConstAtomIterator localIter;
      j = 0;
      for (begin(localIter); localIter.notEnd(); ++localIter) {
         ++j;
         ptr = find(localIter->id());
         if (ptr == 0)
            UTIL_THROW("Unable to find local atom returned by iterator"); 
         if (ptr != localIter.get())
            UTIL_THROW("Inconsistent find(localIter->id()"); 
      }
      if (j != nAtom())
         UTIL_THROW("Number from localIterator != nAtom()"); 

      // Count ghost atoms
      ConstGhostIterator ghostIter;
      j = 0;
      for (begin(ghostIter); ghostIter.notEnd(); ++ghostIter) {
         ++j;
         ptr = find(ghostIter->id());
         if (ptr == 0)
            UTIL_THROW("find(ghostIter->id() == 0"); 
      }
      if (j != nGhost())
         UTIL_THROW("Number from ghostIterator != nGhost()"); 

      return true;
   }

   #ifdef UTIL_MPI
   /**
   * Compute, store and return total number of atoms on all processors.
   */
   void AtomStorage::computeNAtomTotal(MPI::Intracomm& communicator)
   {
      int nAtomLocal = nAtom();
      communicator.Reduce(&nAtomLocal, &nAtomTotal_, 1, 
                          MPI::INT, MPI::SUM, 0);
      if (communicator.Get_rank() !=0) {
         nAtomTotal_ = -1;
      }
   }
   #endif

}
#endif
