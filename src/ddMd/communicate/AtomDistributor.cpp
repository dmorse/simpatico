#ifndef DDMD_ATOM_DISTRIBUTOR_CPP
#define DDMD_ATOM_DISTRIBUTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomDistributor.h"
#include "Domain.h"
#include "Buffer.h"
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/ConstAtomIterator.h>

#include <algorithm>

//#define DDMD_ATOM_DISTRIBUTOR_DEBUG

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   AtomDistributor::AtomDistributor() :
      #ifdef UTIL_MPI 
      sendArrays_(),
      sendSizes_(),
      #endif
      boundaryPtr_(0),
      domainPtr_(0),
      storagePtr_(0),
      #ifdef UTIL_MPI
      bufferPtr_(0),
      #endif
      newPtr_(0),
      cacheCapacity_(1024),
      sendCapacity_(0),
      rankMaxSendSize_(0),
      nCachedTotal_(0),
      nSentTotal_(0),
      isAllocated_(false)
   {  setClassName("AtomDistributor"); }

   /*
   * Destructor.
   */
   AtomDistributor::~AtomDistributor() 
   {}

   /*
   * Retain pointers to associated objects (call on all domain processors).
   */
   void AtomDistributor::associate(Domain& domain, Boundary& boundary, 
                                   AtomStorage& storage, Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      storagePtr_ = &storage;
      #ifdef UTIL_MPI
      bufferPtr_ = &buffer;
      #endif
   }

   /*
   * Set cache capacity.
   */
   void AtomDistributor::setCapacity(int cacheCapacity)
   {  
      if (cache_.capacity() > 0) { 
         UTIL_THROW("Attempt to set cacheCapacity after allocation");
      } 
      cacheCapacity_ = cacheCapacity; 
   }

   /*
   * Read cacheCapacity.
   */
   void AtomDistributor::readParameters(std::istream& in)
   {
      // Read parameter file block
      read<int>(in, "cacheCapacity", cacheCapacity_);
   }

   /*
   * Allocate and initialize memory on master processor (private).
   */
   void AtomDistributor::allocate()
   {
      // Preconditions
      if (isAllocated_) {
         UTIL_THROW("Attempt to re-allocate AtomDistributor");
      }
      if (domainPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
      if (!domainPtr_->isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }
      #ifdef UTIL_MPI
      if (!bufferPtr_->isInitialized()) {
         UTIL_THROW("Buffer is not initialized");
      }
      #endif

      int gridSize  = domainPtr_->grid().size();
      #ifdef UTIL_MPI
      sendCapacity_ = bufferPtr_->atomCapacity();
      #endif

      // Default cacheCapacity_ = (# processors)*(max atoms per send)
      if (cacheCapacity_ <= 0) {
         cacheCapacity_ = gridSize * sendCapacity_;
      }

      // Allocate atom cache and reservoir on the master processor. 
      cache_.allocate(cacheCapacity_);
      reservoir_.allocate(cacheCapacity_);

      // Push all atoms onto the reservoir stack, in reverse order.
      for (int i = cacheCapacity_ - 1; i >= 0; --i) {
         reservoir_.push(cache_[i]);
      }

      #ifdef UTIL_MPI
      // Allocate memory for sendArrays_ matrix, and nullify all elements.
      sendArrays_.allocate(gridSize, sendCapacity_);
      sendSizes_.allocate(gridSize);
      for (int i = 0; i < gridSize; ++i) {
         sendSizes_[i] = 0; 
         for (int j = 0; j < sendCapacity_; ++j) {
            sendArrays_(i, j) = 0;      
         }
      }
      #endif

      // Mark this as allocated.
      isAllocated_ = true;
   }

   #ifdef UTIL_MPI
   /*
   * Initialize the send buffer.
   */ 
   void AtomDistributor::setup() 
   {  
      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("AtomDistributor not initialized");
      }
      if (!domainPtr_->isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }
      if (domainPtr_->gridRank() != 0) {
         UTIL_THROW("This is not the master processor");
      }
      #ifdef UTIL_MPI
      if (!bufferPtr_->isInitialized()) {
         UTIL_THROW("Buffer is not initialized");
      }
      #endif

      // Allocate if needed (on first call).
      // Note: isAllocated_ is set true by the allocate() function.
      if (!isAllocated_) {
         allocate();
      }

      // Check post-allocation conditions 
      if (reservoir_.size() != reservoir_.capacity()) {
         UTIL_THROW("Atom reservoir not full in setup");
      }
      int gridSize  = domainPtr_->grid().size();
      for (int i = 0; i < gridSize; ++i) {
         if (sendSizes_[i] != 0) {
            UTIL_THROW("A sendArray size is not zero in setup"); 
         }
      }

      // Clear buffer and counters
      bufferPtr_->clearSendBuffer();
      bufferPtr_->beginSendBlock(Buffer::ATOM);
      nCachedTotal_ = 0;
      nSentTotal_ = 0;
   }
   #endif

   /*
   * Returns address for a new local Atom.
   */ 
   Atom* AtomDistributor::newAtomPtr()
   {
      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
      #ifdef UTIL_MPI
      if (bufferPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
      #endif
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("AtomDistributor is not allocated");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      if (domainPtr_->gridRank() != 0) {
         UTIL_THROW("This is not the master processor");
      }
      if (newPtr_ != 0) {
         UTIL_THROW("A newPtr_ is still active");
      }

      #ifdef UTIL_MPI
      // If reservoir is empty, send buffer to processor with maximum sendSize_
      if (reservoir_.size() == 0) {

         Atom* ptr;
         int rank  = rankMaxSendSize_;
         int size  = sendSizes_[rank];
         int nSend = std::min(size, sendCapacity_);
         int begin = size - nSend;

         // Pack atoms into buffer, and return pointers for reuse.
         for (int i = begin; i < size; ++i) {
            ptr = sendArrays_(rank, i);
            ptr->packAtom(*bufferPtr_);

            // Return pointer to the reservoir and remove it from sendArrays_.
            reservoir_.push(*sendArrays_(rank, i));
            sendArrays_(rank, i) = 0;
         }
         bool isComplete = false;
         bufferPtr_->endSendBlock(isComplete);
         nSentTotal_ += sendSizes_[rank];
         sendSizes_[rank] = begin;

         // Send the buffer
         bufferPtr_->send(domainPtr_->communicator(), rank);

         // Reinitialize the buffer for reuse.
         bufferPtr_->clearSendBuffer();
         bufferPtr_->beginSendBlock(Buffer::ATOM);
      }
      #endif

      if (reservoir_.size() == 0) {
         UTIL_THROW("Empty cache reservoir (This should not happen)");
      }
 
      // Pop pointer to new atom from reservoir and return that pointer.
      newPtr_ = &reservoir_.pop();
      newPtr_->clear();
      return newPtr_;
   }

   /*
   * Add an atom to the list to be sent.
   */ 
   int AtomDistributor::addAtom() 
   {
      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
      if (boundaryPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("AtomDistributor is not allocated");
      }
      if (domainPtr_->gridRank() != 0) {
         UTIL_THROW("This is not the master processor");
      }
      if (storagePtr_->isCartesian()) {
         UTIL_THROW("Atom coordinates are Cartesian");
      }
      if (newPtr_ == 0) {
         UTIL_THROW("No active newPtr_");
      }

      // Shift position to lie within primary unit cell.
      boundaryPtr_->shiftGen(newPtr_->position());

      #ifdef UTIL_MPI
      // Identify rank of processor that owns this atom.
      int rank = domainPtr_->ownerRank(newPtr_->position());
      #else
      int rank = 0;
      #endif

      // If owned by the master, add this atom to storage.
      // If not owned by the master, queue this atom for sending.
      if (rank == 0) {

         Atom* ptr = storagePtr_->newAtomPtr();
         *ptr = *newPtr_;
         storagePtr_->addNewAtom();
         reservoir_.push(*newPtr_);

         // Note: Atom is returned to reservoir in a dirty state.
         // Atoms must thus be cleared when popped from reservoir.
      }
      #ifdef UTIL_MPI
      else { // if rank !=0

         // Add newPtr_ to array of pointers for processor rank.
         assert(sendSizes_[rank] < sendCapacity_);
         sendArrays_(rank, sendSizes_[rank]) = newPtr_;
         ++sendSizes_[rank];
         ++nCachedTotal_;

         // Check if sendSize for this rank is the largest.
         if (rank != rankMaxSendSize_) {
            if (sendSizes_[rank] > sendSizes_[rankMaxSendSize_]) {
               rankMaxSendSize_ = rank;
            }
         }
 
         // If buffer for the relevant processor is full, send it now.
         if (sendSizes_[rank] == sendCapacity_) {

            // Pack atoms into buffer, and return pointers for reuse.
            Atom* ptr;
            for (int i = 0; i < sendCapacity_; ++i) {
               ptr = sendArrays_(rank, i);
               ptr->packAtom(*bufferPtr_);

               // Push pointer onto reservoir and remove it from sendArrays_.
               reservoir_.push(*sendArrays_(rank, i));
               sendArrays_(rank, i) = 0;
            }
            bool isComplete = false;
            bufferPtr_->endSendBlock(isComplete);
            nSentTotal_ += sendCapacity_;
            sendSizes_[rank] = 0;

            // Send the buffer
            bufferPtr_->send(domainPtr_->communicator(), rank);

            // Reinitialize the send buffer for reuse.
            bufferPtr_->clearSendBuffer();
            bufferPtr_->beginSendBlock(Buffer::ATOM);
         }

      }
      #endif

      // Nullify newPtr_ to release for reuse.
      newPtr_ = 0;

      #ifdef DDMD_ATOM_DISTRIBUTOR_DEBUG
      // Check allocation of cache atoms
      int sendSizeSum = 0;
      int gridSize  = domainPtr_->grid().size();
      for (int i = 0; i < gridSize; ++i) {
         sendSizeSum += sendSizes_[i]; 
      }
      if (sendSizeSum + reservoir_.size() != cacheCapacity_) {
         UTIL_THROW("Error: Inconsistent cache atom count");
      }
      #endif

      // Return rank of processor that owns this atom.
      return rank;
   }

   #ifdef UTIL_MPI
   /*
   * Send any atoms that have not be sent previously.
   *
   * Called only on the master processor.
   */
   void AtomDistributor::send()
   {

      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
      if (boundaryPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      if (domainPtr_->gridRank() != 0) {
         UTIL_THROW("This is not the master processor");
      }
      if (newPtr_ != 0) {
         UTIL_THROW("A newPtr_ is still active");
      }

      int gridSize = domainPtr_->grid().size();
      int i, j;
      bool isComplete = true;
      for (i = 1; i < gridSize; ++i) {

         // Pack all remaining atoms for this processor
         Atom* ptr;
         for (j = 0; j < sendSizes_[i]; ++j) {
            ptr = sendArrays_(i, j);
            ptr->packAtom(*bufferPtr_);

            // Return pointer to atom to the reservoir.
            reservoir_.push(*sendArrays_(i, j));
         }
         bufferPtr_->endSendBlock(isComplete);
         nSentTotal_ += sendSizes_[i];

         // Send buffer, with isComplete flag as true.
         bufferPtr_->send(domainPtr_->communicator(), i);

         // Clear buffer and initialize for resending.
         bufferPtr_->clearSendBuffer();
         if (i < gridSize - 1) {
            bufferPtr_->beginSendBlock(Buffer::ATOM);
         }

         // Reset send arrays to empty state.
         for (int j = 0; j < sendCapacity_; ++j) {
            sendArrays_(i, j) = 0;
         }
         sendSizes_[i] = 0;
      }

      // Compute total number of atoms on all processors.
      // Note: Matching call at end of AtomDistributor::receive()
      storagePtr_->unsetNAtomTotal();
      storagePtr_->computeNAtomTotal(domainPtr_->communicator());

      // Postconditions
      if (reservoir_.size() != reservoir_.capacity()) {
         UTIL_THROW("atomReservoir not empty after final send");
      }
      if (nCachedTotal_ != nSentTotal_) {
         UTIL_THROW("Number cached atoms != number sent");
      }
      if (storagePtr_->nAtomTotal() != nSentTotal_ + storagePtr_->nAtom()) {
         UTIL_THROW("Number atoms received != number sent + nAtom on master");
      }

   }
   #endif

   #ifdef UTIL_MPI
   /*
   * Receive all atoms sent by the master processor.
   *
   * Called by all domain processors except the master.
   */ 
   void AtomDistributor::receive()
   {
      Atom* ptr;                 // Ptr to atom for storage
      int   rank;                // Rank of this processor
      const int source = 0;      // Rank of source processor
      bool  isComplete = false;  // Have all atoms been received?

      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      rank = domainPtr_->gridRank();
      if (rank == 0) {
         UTIL_THROW("AtomDistributor::receive() called on master processor");
      }

      while (!isComplete) {

         // Receive a buffer
         bufferPtr_->recv(domainPtr_->communicator(), source);

         // Unpack the buffer
         isComplete = bufferPtr_->beginRecvBlock();
         while (bufferPtr_->recvSize() > 0) {
            ptr = storagePtr_->newAtomPtr();
            ptr->unpackAtom(*bufferPtr_);
            storagePtr_->addNewAtom();
            if (domainPtr_->ownerRank(ptr->position()) != rank) {
               UTIL_THROW("Error: Atom on wrong processor");
            }
         }
         bufferPtr_->endRecvBlock();

      }

      // Compute total number of atoms on all processors.
      // Note: Matching call at end of AtomDistributor::send()
      storagePtr_->unsetNAtomTotal();
      storagePtr_->computeNAtomTotal(domainPtr_->communicator());
   }
   #endif

   /*
   * Validate distribution of atoms, return total number of atoms.
   * Called on all processors. Correct return value only on master.
   */
   int AtomDistributor::validate()
   {
      storagePtr_->isValid(domainPtr_->communicator());

      // Check that number of atoms = nSentTotal
      int nAtomTotal = 0;
      storagePtr_->computeNAtomTotal(domainPtr_->communicator());
      if (domainPtr_->isMaster()) {
         nAtomTotal = storagePtr_->nAtomTotal();
         if (nAtomTotal != nSentTotal_ + storagePtr_->nAtom()) {
            UTIL_THROW("nAtomTotal != nAtom after distribution");
         }
      }

      // Check that every atom is in correct processor domain.
      // Coordinates must scaled / generalized, rather than Cartesian.
      double coordinate;
      ConstAtomIterator atomIter;
      int i;
      storagePtr_->begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         for (i=0; i < Dimension; ++i) {
            coordinate = atomIter->position()[i];
            if (coordinate < domainPtr_->domainBound(i, 0)) {
               UTIL_THROW("coordinate < domainBound(i, 0)");
            }
            if (coordinate >= domainPtr_->domainBound(i, 1)) {
               UTIL_THROW("coordinate >= domainBound(i, 1)");
            }
         }
      } 

      return nAtomTotal; // Only valid on master, returns 0 otherwise.
   }

} // namespace DdMd
#endif
