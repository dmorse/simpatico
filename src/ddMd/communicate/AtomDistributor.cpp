#ifndef ATOM_DISTRIBUTOR_CPP
#define ATOM_DISTRIBUTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomDistributor.h"
#include "Domain.h"
#include "Buffer.h"
#include <ddMd/storage/AtomStorage.h>

#include <algorithm>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   AtomDistributor::AtomDistributor() 
    : sendArrays_(),
      sendSizes_(),
      boundaryPtr_(0),
      domainPtr_(0),
      bufferPtr_(0),
      newPtr_(0),
      cacheCapacity_(0),
      sendCapacity_(0),
      rankMaxSendSize_(0),
      nCachedTotal_(0),
      nSentTotal_(0)
   {}

   /*
   * Destructor.
   */
   AtomDistributor::~AtomDistributor() 
   {}

   /*
   * Retain pointers to associated objects.
   */
   void AtomDistributor::associate(Domain& domain, Boundary& boundary, 
                                   Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      bufferPtr_ = &buffer;
   }

   /*
   * Set cache capacity and allocate all required memory.
   */
   void AtomDistributor::setParam(int cacheCapacity)
   {
      cacheCapacity_ = cacheCapacity;
      allocate();
   }

   /*
   * Read cacheCapacity and allocate all required memory.
   */
   void AtomDistributor::readParam(std::istream& in)
   {
      // Read parameter file block
      readBegin(in, "AtomDistributor");
      read<int>(in, "cacheCapacity", cacheCapacity_);
      readEnd(in);
 
      // Do actual allocation
      allocate();
   }

   /*
   * Allocate memory and initialize state (private method).
   */
   void AtomDistributor::allocate()
   {
      // Preconditions
      if (bufferPtr_ == 0) {
         UTIL_THROW("AtomDistributor not initialized");
      }
      if (!bufferPtr_->isInitialized()) {
         UTIL_THROW("Buffer not initialized");
      }
      if (!domainPtr_->isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }

      int gridSize  = domainPtr_->grid().size();
      int rank      = domainPtr_->gridRank();
      sendCapacity_ = bufferPtr_->atomCapacity();

      // If master processor
      if (rank == 0) {
         int i, j;

         // Default cacheCapacity_ = (# processors)*(max atoms per send)
         if (cacheCapacity_ <= 0) {
            cacheCapacity_ = gridSize * sendCapacity_;
         }

         // Allocate memory for array of atoms on the master processor. 
         cache_.allocate(cacheCapacity_);
         reservoir_.allocate(cacheCapacity_);

         // Push all atoms onto the reservoir stack, in reverse order.
         for (i = cacheCapacity_ - 1; i >= 0; --i) {
            reservoir_.push(cache_[i]);
         }

         // Allocate memory for sendArrays_ matrix, and nullify all elements.
         sendArrays_.allocate(gridSize, sendCapacity_);
         sendSizes_.allocate(gridSize);
         for (i = 0; i < gridSize; ++i) {
            sendSizes_[i] = 0; 
            for (j = 0; j < sendCapacity_; ++j) {
               sendArrays_(i, j) = 0;      
            }
         }

      }

   }

   #ifdef UTIL_MPI

   /*
   * Initialize the send buffer.
   */ 
   void AtomDistributor::initSendBuffer() 
   {  
      bufferPtr_->clearSendBuffer(); 
      bufferPtr_->beginSendBlock(Buffer::ATOM); 
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
      if (bufferPtr_ == 0) {
         UTIL_THROW("AtomDistributor is not initialized");
      }
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

         int rank  = rankMaxSendSize_;
         int size  = sendSizes_[rank];
         int nSend = std::min(size, sendCapacity_);
         int begin = size - nSend;

         // Pack atoms into buffer, and return pointers for reuse.
         for (int i = begin; i < size; ++i) {
            bufferPtr_->packAtom(*sendArrays_(rank, i));

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

      // Return pointer to new atom.
      newPtr_ = &reservoir_.pop();
      return newPtr_;

   }

   /*
   * Add an atom to the list to be sent.
   */ 
   int AtomDistributor::addAtom(AtomStorage& storage) 
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
      if (newPtr_ == 0) {
         UTIL_THROW("No active newPtr_");
      }

      // Shift position to lie within boundary.
      boundaryPtr_->shift(newPtr_->position());

      // Identify rank of processor that owns this atom.
      int rank = domainPtr_->ownerRank(newPtr_->position());

      // If not owned by the master, queue this atom for sending.
      if (rank != 0) {

         #ifdef UTIL_MPI

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
            for (int i = 0; i < sendCapacity_; ++i) {
               bufferPtr_->packAtom(*sendArrays_(rank, i));

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

         #else 

         UTIL_THROW("Atom not owned by master but UTIL_MPI not defined");

         #endif

      } else { // rank == 0

         Atom* ptr  = storage.newAtomPtr();
         *ptr = *newPtr_;
         storage.addNewAtom();

         reservoir_.push(*newPtr_); 
      }

      // Nullify newPtr_ to release for reuse.
      newPtr_ = 0;

      // Return rank of processor that owns this atom.
      return rank;

   }

   /*
   * Send any atoms that have not be sent previously.
   *
   * This method should be called only by the master processor.
   */
   void AtomDistributor::send()
   {

      #ifdef UTIL_MPI

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
         for (j = 0; j < sendSizes_[i]; ++j) {
            bufferPtr_->packAtom(*sendArrays_(i, j));

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

         // Reset atomPtrs array to empty state.
         for (int j = 0; j < sendCapacity_; ++j) {
            sendArrays_(i, j) = 0;
         }
         sendSizes_[i] = 0;
      }

      // Postconditions
      if (reservoir_.size() != reservoir_.capacity()) {
         UTIL_THROW("atomReservoir not empty after final send");
      }
      if (nCachedTotal_ != nSentTotal_) {
         UTIL_THROW("Number of cached atoms != number sent");
      }

      nCachedTotal_ = 0;
      nSentTotal_   = 0;
      #endif

   }

   /*
   * Receive all atoms sent by the master processor.
   *
   * Called by all processors except the master.
   */ 
   void AtomDistributor::receive(AtomStorage& storage)
   {
      #ifdef UTIL_MPI
      Atom* ptr;                 // Ptr to atom for storage
      const int source = 0;      // Rank of source processor
      int   rank;                // Rank of this processor
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
            ptr = storage.newAtomPtr();
            bufferPtr_->unpackAtom(*ptr);
            storage.addNewAtom();
            if (domainPtr_->ownerRank(ptr->position()) != rank) {
               UTIL_THROW("Error: Atom on wrong processor");
            }
         }

      }
      #endif
   }

}
#endif
