#ifndef DDMD_GROUP_COLLECTOR_TPP
#define DDMD_GROUP_COLLECTOR_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupCollector.h"
#include "Domain.h"
#include "Buffer.h"
#include <ddMd/storage/GroupStorage.tpp>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int N>
   GroupCollector<N>::GroupCollector() 
    : domainPtr_(0),
      storagePtr_(0),
      bufferPtr_(0),
      source_(-1),
      recvArrayCapacity_(256),
      recvBufferSize_(-1),
      recvArraySize_(-1),
      recvArrayId_(-1),
      isComplete_(false)
   {  setClassName("GroupCollector"); }

   /*
   * Destructor.
   */
   template <int N>
   GroupCollector<N>::~GroupCollector() 
   {}

   /*
   * Retain pointers to associated objects.
   */
   template <int N>
   void GroupCollector<N>::associate(Domain& domain, GroupStorage<N>& storage, 
                                 Buffer& buffer)
   {
      domainPtr_  = &domain;
      storagePtr_ = &storage;
      bufferPtr_  = &buffer;
   }

   /*
   * Set recvArray capacity (only needed on master).
   */
   template <int N>
   void GroupCollector<N>::setCapacity(int recvArrayCapacity)
   {  
      if (recvArrayCapacity <= 0) {
         UTIL_THROW("Attempt to set nonpositive recvArrayCapacity");
      }  
      if (recvArray_.capacity() > 0) { 
         UTIL_THROW("Attempt to set recvArrayCapacity after allocation");
      } 
      recvArrayCapacity_ = recvArrayCapacity; 
   }

   /*
   * Setup on master processor just before loop over groups.
   */
   template <int N>
   void GroupCollector<N>::setup()
   {
      // Preconditions
      if (!domainPtr_) {
         UTIL_THROW("Collector not initialized: No associated domain");
      }
      if (!domainPtr_->isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }
      if (!domainPtr_->isMaster()) {
         UTIL_THROW("Not the master processor");
      }
      if (!bufferPtr_->isInitialized()) {
         UTIL_THROW("Buffer not allocated");
      }

      // Allocate recvArray if necessary
      if (recvArray_.capacity() == 0) {
         if (recvArrayCapacity_ == 0) {
             UTIL_THROW("recvArrayCapacity_ not set");
         }
         recvArray_.allocate(recvArrayCapacity_); 
      }

      source_ = 0;          // rank of source node
      recvBufferSize_ = 0;  // number of groups in MPI buffer 
      recvArraySize_ = 0;   // number of groups in recvArray_ 
      recvArrayId_ = 0;     // id of current group in recvArray_
      isComplete_ = false;  // not finished with current processor

      // Initialize Group iterator on master processor.
      storagePtr_->begin(iterator_);
   }

   #ifdef UTIL_MPI
   /*
   * Returns address for a new Group.
   */ 
   template <int N>
   Group<N>* GroupCollector<N>::nextPtr()
   {
      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("GroupCollector has not been initialized");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      if (!domainPtr_->isMaster()) {
         UTIL_THROW("Not the master processor");
      }
      if (recvArray_.capacity() <= 0) {
         UTIL_THROW("Cache not allocated");
      }

      // If master processor
      Group<N>* groupPtr;
      Atom*     atomPtr;
      if (source_ == 0) {
         while (!isComplete_) {
            if (iterator_.notEnd()) {
               groupPtr = iterator_.get();
               ++iterator_;
               atomPtr = groupPtr->atomPtr(0);
               if (atomPtr) {
                  if (!atomPtr->isGhost()) {
                     return groupPtr;
                  }
               }
            } else {
               recvBufferSize_ = 0;
               recvArraySize_ = 0;
               recvArrayId_ = 0;
               isComplete_ = true;
            }
         }
      }
     
      // While at end of recvArray_, or while array is empty.
      while (recvArrayId_ == recvArraySize_) {

         // If receive buffer is empty
         if (recvBufferSize_ == 0) {

            // If processing of items from processor source_ is complete.
            if (isComplete_) {
               ++source_;
               recvArraySize_ = 0;
               recvArrayId_ = 0;
               isComplete_ = false;
               // If last processor is complete, return null pointer.
               if (source_ == domainPtr_->grid().size()) {
                  source_ = 0;
                  return 0;
               }
            }

            // Send request to processor source_ 
            int message = source_;
            domainPtr_->communicator().Send(&message, 1, MPI::INT, 
                                               source_, message);
      
            // Receive buffer from processor source_
            bufferPtr_->recv(domainPtr_->communicator(), source_);
            isComplete_ = bufferPtr_->beginRecvBlock();
            recvBufferSize_ = bufferPtr_->recvSize();
         }

         // Unpack groups from buffer into recvArray_.
         if (recvBufferSize_ > 0) {
            recvArraySize_ = 0;
            recvArrayId_ = 0;
            while (bufferPtr_->recvSize() > 0 
                   && recvArraySize_ < recvArray_.capacity()) 
            {
               recvArray_[recvArraySize_].unpack(*bufferPtr_);
               ++recvArraySize_;
               --recvBufferSize_;
               if (recvBufferSize_ != bufferPtr_->recvSize()) {
                  UTIL_THROW("Inconsistent buffer receive counters");
               }
            }
         }

      }

      // Return current item from recvArray.
      ++recvArrayId_;
      return &recvArray_[recvArrayId_ - 1];

   }

   /*
   * Send all groups from this process.
   *
   * Call on every processor except the master.
   */
   template <int N>
   void GroupCollector<N>::send()
   {

      // Preconditions
      if (!domainPtr_) {
         UTIL_THROW("Collector not initialized: null domainPtr_");
      }
      if (!domainPtr_->isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }
      if (domainPtr_->isMaster()) {
         UTIL_THROW("GroupCollector<N>::send() called from master node.");
      }
      if (!storagePtr_) {
         UTIL_THROW("Collector not initialized: Null storagePtr_");
      }
      if (storagePtr_->capacity() <= 0) {
         UTIL_THROW("GroupStorage not initialized");
      }
      if (!bufferPtr_) {
         UTIL_THROW("Collector not initialized: Null bufferPtr_");
      }
      if (!bufferPtr_->isInitialized()) {
         UTIL_THROW("Buffer not allocated");
      }

      Atom* atomPtr = 0;
      int message;
      int tag;
      int bufferCapacity = bufferPtr_->groupCapacity<N>();

      // Initialize group iterator
      storagePtr_->begin(iterator_);

      isComplete_= false;
      while (!isComplete_) {

         // Receive notice from master to send groups (blocking receive)
         tag = domainPtr_->communicator().Get_rank();
         domainPtr_->communicator().Recv(&message, 1, MPI::INT, 0, tag);

         // Pack buffer with groups
         int recvArraySize_ = 0;
         isComplete_ = iterator_.isEnd();
         bufferPtr_->clearSendBuffer();
         bufferPtr_->beginSendBlock(Buffer::GROUP2 +  N - 2);
         while (recvArraySize_ < bufferCapacity && !isComplete_) {
            // Get pointer to first atom in Group
            // Send group only if this is a local atom.
            atomPtr = iterator_->atomPtr(0);
            if (atomPtr) {
               if (!atomPtr->isGhost()) {
                  iterator_->pack(*bufferPtr_);
                  ++recvArraySize_;
               }
            }
            ++iterator_;
            isComplete_ = iterator_.isEnd();
         }
         bufferPtr_->endSendBlock(isComplete_);

         // Send buffer to master
         bufferPtr_->send(domainPtr_->communicator(), 0);

      }

   }
   #endif

}
#endif // ifndef DDMD_GROUP_COLLECTOR_TPP
