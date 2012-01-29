#ifndef COLLECTOR_CPP
#define COLLECTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Collector.h"
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
   Collector::Collector() 
    : domainPtr_(0),
      bufferPtr_(0),
      source_(-1),
      recvArraySize_(-1),
      recvArrayId_(-1),
      isComplete_(false)
   {}

   /*
   * Destructor.
   */
   Collector::~Collector() 
   {}

   /*
   * Retain pointers to associated objects.
   */
   void Collector::initialize(AtomStorage& storage, Domain& domain, 
                              Buffer& buffer)
   {
      // Preconditions
      if (!buffer.isInitialized()) {
         UTIL_THROW("Buffer not allocated");
      }
      if (!domain.isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }
      if (!domain.isMaster()) {
         UTIL_THROW("Not the master processor");
      }

      // Retain pointers to domain and buffer
      domainPtr_  = &domain;
      bufferPtr_  = &buffer;

      // Allocate memory for array of atoms on the master processor. 
      recvArray_.allocate(buffer.atomCapacity());

      source_     = 0;      // rank of source node
      recvArraySize_      = 0;      // number of atoms in recvArray_ array
      recvArrayId_     = 0;      // id of current atom in recvArray_
      isComplete_ = false;

      // Initialize Atom iterator on master processor.
      storage.begin(iterator_);

   }

   #ifdef UTIL_MPI
   /*
   * Returns address for a new Atom.
   */ 
   Atom* Collector::nextPtr()
   {
      // Preconditions
      if (domainPtr_ == 0) {
         UTIL_THROW("Collector has not been initialized");
      }
      if (!domainPtr_->isInitialized() != 0) {
         UTIL_THROW("Domain is not initialized");
      }
      if (!domainPtr_->isMaster()) {
         UTIL_THROW("Not the master processor");
      }

      // If master processor
      Atom* ptr;
      if (source_ == 0) {
         if (!iterator_.atEnd()) {
            ptr = iterator_.get();
            ++iterator_;
            return ptr;
         } else {
            recvArrayId_ = 0;
            recvArraySize_ = 0;
            isComplete_ = true;
         }
      }
     
      // If at end of recvArray_ array, or array is empty.
      if (recvArrayId_ == recvArraySize_) {

         // If processing of items from processor source_ is complete.
         if (isComplete_) {
            ++source_;
            // If last processor is complete, return null pointer.
            if (source_ == domainPtr_->grid().size()) {
               return 0;
            } else {
               recvArrayId_ = 0;
               recvArraySize_  = 0;
               isComplete_ = false;
            }
         }

         // Send message to processor source_ requesting a buffer
         int message = source_;
         domainPtr_->communicator().Send(&message, 1, MPI::INT, source_, message);

         // Receive buffer from processor source_
         bufferPtr_->recv(domainPtr_->communicator(), source_);

         // Unpack all atoms from buffer into recvArray_.
         recvArraySize_ = 0;
         isComplete_ = bufferPtr_->beginRecvBlock();
         while (bufferPtr_->recvSize() > 0) {
            bufferPtr_->unpackAtom(recvArray_[recvArraySize_]);
            ++recvArraySize_;
         }

         recvArrayId_ = 0;

      }

      // Return current item from recvArray.
      ++recvArrayId_;
      return &recvArray_[recvArrayId_ - 1];


   }

   /*
   * Send all atoms from this process.
   *
   * Call on every processor except the master.
   */
   void Collector::send(AtomStorage& storage, Domain& domain, Buffer& buffer)
   {

      // Preconditions
      if (!buffer.isInitialized()) {
         UTIL_THROW("Buffer not allocated");
      }
      if (!domain.isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }
      if (domain.isMaster()) {
         UTIL_THROW("Collector::send() called from master node.");
      }

      // Initialize atom iterator
      storage.begin(iterator_);

      isComplete_= false;
      while (!isComplete_) {

         // Receive notice from master to send atoms (blocking receive)
         int message;
         int tag = domain.communicator().Get_rank();
         domain.communicator().Recv(&message, 1, MPI::INT, 0, tag);

         // Pack buffer with atoms
         int recvArraySize_ = 0;
         isComplete_ = iterator_.atEnd();
         buffer.clearSendBuffer();
         buffer.beginSendBlock(Buffer::ATOM);
         while (recvArraySize_ < buffer.atomCapacity() && !isComplete_) {
            buffer.packAtom(*iterator_);
            ++recvArraySize_;
            ++iterator_;
            isComplete_ = iterator_.atEnd();
         }
         buffer.endSendBlock(isComplete_);

         // Send buffer to master
         buffer.send(domain.communicator(), 0);

      }

   }
   #endif

}
#endif
