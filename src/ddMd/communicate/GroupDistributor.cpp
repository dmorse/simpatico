#ifndef GROUP_DISTRIBUTOR_CPP
#define GROUP_DISTRIBUTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupDistributor.h"
#include "Domain.h"
#include "Buffer.h"
#include <ddMd/storage/GroupStorage.h>

#include <algorithm>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int N>
   GroupDistributor<N>::GroupDistributor() 
    : sendArrays_(),
      sendSizes_(),
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
   template <int N>
   GroupDistributor<N>::~GroupDistributor() 
   {}

   /*
   * Retain pointers to associated objects.
   */
   template <int N>
   void GroupDistributor<N>::associate(Domain& domain, Buffer& buffer)
   {
      domainPtr_    = &domain;
      bufferPtr_    = &buffer;
   }

   /*
   * Set cache capacity and allocate all required memory.
   */
   template <int N>
   void GroupDistributor<N>::setParam(int cacheCapacity)
   {
      cacheCapacity_ = cacheCapacity;
      allocate();
   }

   /*
   * Read cacheCapacity and allocate all required memory.
   */
   template <int N>
   void GroupDistributor<N>::readParam(std::istream& in)
   {
      // Read parameter file block
      readBegin(in, "GroupDistributor");
      read<int>(in, "cacheCapacity", cacheCapacity_);
      readEnd(in);
 
      // Do actual allocation
      allocate();
   }

   /*
   * Allocate memory and initialize state (private method).
   */
   template <int N>
   void GroupDistributor<N>::allocate()
   {
      // Preconditions
      if (bufferPtr_ == 0) {
         UTIL_THROW("GroupDistributor not initialized");
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
   template <int N>
   void GroupDistributor<N>::initSendBuffer() 
   {  
      bufferPtr_->clearSendBuffer(); 
      bufferPtr_->beginSendBlock(Buffer::ATOM); 
   }

   #endif

   /*
   * Returns address for a new local Group.
   */ 
   template <int N>
   Group<N>* GroupDistributor<N>::newPtr()
   {
      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (bufferPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("GroupDistributor is not allocated");
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

         // Pack groups into buffer, and return pointers for reuse.
         for (int i = begin; i < size; ++i) {
            bufferPtr_->packGroup(*sendArrays_(rank, i));

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
   template <int N>
   void 
   GroupDistributor<N>::add(GroupStorage<N>& storage, DArray<int> atomOwners)
   {
      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("GroupDistributor is not allocated");
      }
      if (domainPtr_->gridRank() != 0) {
         UTIL_THROW("This is not the master processor");
      }
      if (newPtr_ == 0) {
         UTIL_THROW("No active newPtr_");
      }

      int  ranks[N];
      int  nRanks = 0;
      int  i, j, k;
      bool ownedByMaster = false;
      bool newOwner;
      for (i = 0; i < N; ++i) {
         k = newPtr_->atomId(i);
         if (k = 0) {
            ownedByMaster = true;
         } else {
            newOwner = true;
            for (j = 0; j < nRanks; ++j) {
               if (k = ranks[j]) {
                  newOwner = false;
               }
            }
            if (newOwner) {
               ranks[nRanks] = k;
               ++nRanks;
            }
         }
      }

      // If not owned by the master, queue this atom for sending.
      if (ownedByMaster) {

         Group<N>* ptr  = storage.newPtr();
         *ptr = *newPtr_;
         storage.addNewGroup();

         reservoir_.push(*newPtr_); 

      #ifndef UTIL_MPI
      } else {
         UTIL_THROW("Group not owned by master but UTIL_MPI not defined");
      #endif
      }
    
      #ifdef UTIL_MPI
      for (i = 0; i < nRanks; ++i) {

         rank = ranks[i];

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
               bufferPtr_->packGroup(*sendArrays_(rank, i));

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

      // Nullify newPtr_ to release for reuse.
      newPtr_ = 0;


      // Return rank of processor that owns this atom.
      return rank;
      #endif

   }

   /*
   * Send any atoms that have not be sent previously.
   *
   * This method should be called only by the master processor.
   */
   template <int N>
   void GroupDistributor<N>::send()
   {

      #ifdef UTIL_MPI

      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
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
            bufferPtr_->packGroup(*sendArrays_(i, j));

            // Return pointer to atom to the reservoir.
            reservoir_.push(*sendArrays_(i, j));

         }
         bufferPtr_->endSendBlock(isComplete);
         nSentTotal_ += sendSizes_[i];

         // Send buffer, with isComplete flag as true.
         bufferPtr_->send(domainPtr_->communicator(), i);

         // Clear buffer and initialize for resending.
         bufferPtr_->clearSendBuffer();
         bufferPtr_->beginSendBlock(Buffer::ATOM);

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
   template <int N>
   void GroupDistributor<N>::receive(GroupStorage<N>& storage)
   {

      #ifdef UTIL_MPI

      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      if (domainPtr_->gridRank() == 0) {
         UTIL_THROW("GroupDistributor<N>::receive() called on master processor");
      }

      Group<N>* ptr;
      const int source = 0;
      bool  isComplete = false;
      while (!isComplete) {

         // Receive a buffer
         bufferPtr_->recv(domainPtr_->communicator(), source);

         // Unpack the buffer
         isComplete = bufferPtr_->beginRecvBlock();
         while (bufferPtr_->recvSize() > 0) {
            ptr = storage.newPtr();
            bufferPtr_->unpackGroup(*ptr);
            storage.add();
         }

      }

      #endif

   }

}
#endif
