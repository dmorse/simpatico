/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomCollector.h"
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
   AtomCollector::AtomCollector() 
    : domainPtr_(0),
      bufferPtr_(0),
      source_(-1),
      recvBufferSize_(-1),
      recvArrayCapacity_(256),
      recvArraySize_(-1),
      recvArrayId_(-1),
      isComplete_(false)
   {  setClassName("AtomCollector"); }

   /*
   * Destructor.
   */
   AtomCollector::~AtomCollector() 
   {}

   /*
   * Retain pointers to associated objects.
   *
   * Call on all domain notes.
   */
   void AtomCollector::associate(Domain& domain, AtomStorage& storage, 
                                 Buffer& buffer)
   {
      domainPtr_  = &domain;
      storagePtr_ = &storage;
      bufferPtr_  = &buffer;
   }

   /*
   * Set recvArray cache capacity.
   */
   void AtomCollector::setCapacity(int recvArrayCapacity)
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
   * Setup before receiving loop - call only on master.
   */
   void AtomCollector::setup()
   {
      // Preconditions
      if (!domainPtr_) {
         UTIL_THROW("Collector not initialized");
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

      // Allocate recvArray if not done previously
      if (recvArray_.capacity() == 0) {
         if (recvArrayCapacity_ == 0) {
            UTIL_THROW("recvArrayCapacity_ not set");
         }
         recvArray_.allocate(recvArrayCapacity_);
      }

      source_ = 0;          // rank of source node
      recvBufferSize_ = 0;  // number of atoms in MPI buffer 
      recvArraySize_ = 0;   // number of atoms in recvArray_ 
      recvArrayId_ = 0;     // id of current atom in recvArray_
      isComplete_ = false;  // not finished with current processor

      // Initialize Atom iterator on master processor.
      storagePtr_->begin(iterator_);
   }

   #ifdef UTIL_MPI
   /*
   * Return address for a new Atom, or null when all are received.
   *
   * Call this function only on the master processor.
   */ 
   Atom* AtomCollector::nextPtr()
   {
      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("AtomCollector has not been initialized");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      if (!domainPtr_->isMaster()) {
         UTIL_THROW("Not the master processor");
      }

      // If still processing atoms from master processor
      Atom* ptr;
      if (source_ == 0) {
         if (iterator_.notEnd()) {
            ptr = iterator_.get();
            ++iterator_;
            return ptr;
         } else {
            recvBufferSize_ = 0;
            recvArraySize_ = 0;
            recvArrayId_ = 0;
            isComplete_ = true;
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

            assert(recvBufferSize_ == 0); // recv buffer is empty
            assert(!isComplete_);         // source_ is not completed
            
            // Send request to processor source_ .
            int message = source_;
            MPI_Send(&message, 1, MPI_INT, source_, message,
                     domainPtr_->communicator());
      
            // Receive buffer from processor source_
            bufferPtr_->recv(domainPtr_->communicator(), source_);
            isComplete_ = bufferPtr_->beginRecvBlock();
            recvBufferSize_ = bufferPtr_->recvSize();
         }

         // Unpack atoms from recv buffer into recvArray_.
         if (recvBufferSize_ > 0) {
            recvArraySize_ = 0;
            recvArrayId_ = 0;
            if (recvBufferSize_ != bufferPtr_->recvSize()) {
               UTIL_THROW("Inconsistent buffer receive counters");
            }
            while (bufferPtr_->recvSize() > 0 
                   && recvArraySize_ < recvArray_.capacity()) 
            {
               recvArray_[recvArraySize_].unpackAtom(*bufferPtr_);
               ++recvArraySize_;
               --recvBufferSize_;
               if (recvBufferSize_ != bufferPtr_->recvSize()) {
                  UTIL_THROW("Inconsistent buffer receive counters");
               }
            }
            if (bufferPtr_->recvSize() == 0) {
               bufferPtr_->endRecvBlock();
            }
         }

      }

      // Return current item from recvArray.
      ++recvArrayId_;
      return &recvArray_[recvArrayId_ - 1];

   }

   /*
   * Send all atoms from this process.
   *
   * Call on every domain processor except the master.
   */
   void 
   AtomCollector::send()
   {

      // Preconditions
      if (!domainPtr_->isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }
      if (domainPtr_->isMaster()) {
         UTIL_THROW("AtomCollector::send() called from master node.");
      }
      if (!bufferPtr_->isInitialized()) {
         UTIL_THROW("Buffer is not initialized");
      }

      // Initialize atom iterator
      storagePtr_->begin(iterator_);

      isComplete_= false;
      while (!isComplete_) {

         // Receive notice from master to send atoms (blocking receive)
         int tag;
         MPI_Comm_rank(domainPtr_->communicator(), &tag);
         int message;
         MPI_Status status;
         MPI_Recv(&message, 1, MPI_INT, 0, tag, 
                  domainPtr_->communicator(), &status);

         // Pack buffer with atoms
         int recvArraySize_ = 0;
         isComplete_ = iterator_.isEnd();
         bufferPtr_->clearSendBuffer();
         bufferPtr_->beginSendBlock(Buffer::ATOM);
         while (recvArraySize_ < bufferPtr_->atomCapacity() && !isComplete_) {
            iterator_->packAtom(*bufferPtr_);
            ++recvArraySize_;
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
