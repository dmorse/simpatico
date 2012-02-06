#ifndef GROUP_DISTRIBUTOR_CPP
#define GROUP_DISTRIBUTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupDistributor.h"
//#include "Buffer.h"
#include "Domain.h"
#include <ddMd/storage/GroupStorage.h>
#include <ddMd/storage/AtomStorage.h>

#include <algorithm>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int N>
   GroupDistributor<N>::GroupDistributor() 
    : cache_(),
      newPtr_(0),
      atomStoragePtr_(0),
      groupStoragePtr_(0),
      domainPtr_(0),
      bufferPtr_(0),
      cacheCapacity_(0),
      cacheSize_(0),
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
   void GroupDistributor<N>::associate(GroupStorage<N>& groupStorage, 
                                       AtomStorage& atomStorage,
                                       Domain& domain, Buffer& buffer)
   {
      groupStoragePtr_ = &groupStorage;
      atomStoragePtr_ = &atomStorage;
      domainPtr_ = &domain;
      bufferPtr_ = &buffer;
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
      if (atomStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor not initialized");
      }
      cache_.allocate(cacheCapacity_);
   }

   /*
   *
   */
   template <int N>
   void GroupDistributor<N>::initSendBuffer()
   {
      bufferPtr_->clearSendBuffer();
      bufferPtr_->beginSendBlock(Buffer::GROUP, N);
   }

   /*
   * Returns address for a new local Group.
   */ 
   template <int N>
   Group<N>* GroupDistributor<N>::newPtr()
   {
      // Preconditions
      if (atomStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("GroupDistributor is not allocated");
      }
      if (newPtr_ != 0) {
         UTIL_THROW("A newPtr_ is still active");
      }

      // #ifdef UTIL_MPI
      if (cacheSize_ == cacheCapacity_) {
          bool isComplete = false;
          int  source = 0;
          bufferPtr_->endSendBlock(isComplete);
          bufferPtr_->bcast(domainPtr_->communicator(), source);
          nSentTotal_ += cacheSize_;
          cacheSize_ = 0;
          bufferPtr_->clearSendBuffer();
          bufferPtr_->beginSendBlock(Buffer::GROUP, N);
      }
      // #endif
      newPtr_ = &cache_[cacheSize_];
      return newPtr_;
   }

   /*
   * Add an atom to the list to be sent.
   */ 
   template <int N>
   void GroupDistributor<N>::add()
   {
      // Preconditions
      if (atomStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("GroupDistributor is not allocated");
      }

      // If this group has atoms on the master, add to groupStorage.
      if (atomStoragePtr_->groupHasAtoms(*newPtr_)) {

         Group<N>* ptr  = groupStoragePtr_->newPtr();
         *ptr = *newPtr_;
         groupStoragePtr_->add();

      }
    
      // Nullify newPtr_ to release for reuse.
      newPtr_ = 0;
      ++cacheSize_;

   }

   /*
   * Send any atoms that have not be sent previously.
   *
   * This method should be called only by the master processor.
   */
   template <int N>
   void GroupDistributor<N>::send()
   {

      // #ifdef UTIL_MPI
      // Preconditions
      if (atomStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("GroupDistributor is not allocated");
      }
      if (newPtr_ != 0) {
         UTIL_THROW("A newPtr_ is still active");
      }

      if (cacheSize_ > 0) {
          bool isComplete = false;
          int  source = 0;
          bufferPtr_->endSendBlock(isComplete);
          bufferPtr_->bcast(domainPtr_->communicator(), source);
          nSentTotal_ += cacheSize_;
          cacheSize_ = 0;
          bufferPtr_->clearSendBuffer();
          bufferPtr_->beginSendBlock(Buffer::GROUP, N);
      }

      nSentTotal_   = 0;
      // #endif

   }

   /*
   * Receive all atoms sent by the master processor.
   *
   * Called by all processors except the master.
   */ 
   template <int N>
   void GroupDistributor<N>::receive()
   {

      // #ifdef UTIL_MPI
      // Preconditions ??

      Group<N>* ptr;
      const int source = 0;
      bool  isComplete = false;
      bool  myGroup = true;
      while (!isComplete) {

         // Receive broadcast
         bufferPtr_->bcast(domainPtr_->communicator(), source);

         // Unpack the buffer, 
         isComplete = bufferPtr_->beginRecvBlock();
         while (bufferPtr_->recvSize() > 0) {
            if (atomStoragePtr_->groupHasAtoms(*newPtr_)) {
               ptr = groupStoragePtr_->newPtr();
               bufferPtr_->unpackGroup(*ptr);
               groupStoragePtr_->add();
            } else {
               bufferPtr_->discardGroup<N>();
            }
         }

      }
      // #endif

   }

}
#endif
