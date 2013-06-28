#ifndef DDMD_GROUP_DISTRIBUTOR_TPP
#define DDMD_GROUP_DISTRIBUTOR_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupDistributor.h"
#include "Domain.h"
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/GroupStorage.tpp>
#include <ddMd/storage/GroupIterator.h>
#include <ddMd/storage/ConstGroupIterator.h>

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
      domainPtr_(0),
      atomStoragePtr_(0),
      groupStoragePtr_(0),
      bufferPtr_(0),
      sendType_(Buffer::NONE),
      nAtomRecv_(0),
      nSentTotal_(0),
      cacheSize_(0),
      cacheCapacity_(0)
   {  setClassName("GroupDistributor"); }

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
   void GroupDistributor<N>::associate(Domain& domain, 
                                       AtomStorage& atomStorage,
                                       GroupStorage<N>& groupStorage, 
                                       Buffer& buffer)
   {
      domainPtr_ = &domain;
      atomStoragePtr_ = &atomStorage;
      groupStoragePtr_ = &groupStorage;
      bufferPtr_ = &buffer;
   }

   /*
   * Set cache capacity and allocate all required memory.
   */
   template <int N>
   void GroupDistributor<N>::initialize(int cacheCapacity)
   {
      cacheCapacity_ = cacheCapacity;
      allocate();
   }

   /*
   * Read cacheCapacity and allocate all required memory.
   */
   template <int N>
   void GroupDistributor<N>::readParameters(std::istream& in)
   {
      read<int>(in, "cacheCapacity", cacheCapacity_);
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
      nAtomRecv_ = 0;
      newPtr_ = 0;
   }

   /*
   *
   */
   template <int N>
   void GroupDistributor<N>::setup()
   {
      bufferPtr_->clearSendBuffer();
      bufferPtr_->beginSendBlock(Buffer::GROUP2 + N - 2);
      nAtomRecv_ = 0;
      newPtr_ = 0;
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
      if (groupStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (domainPtr_->gridRank() != 0) {
         UTIL_THROW("GroupDistributor::add called on slave node");
      }
      if (cache_.capacity() <= 0) {
         UTIL_THROW("GroupDistributor cache is not allocated");
      }
      if (newPtr_ != 0) {
         UTIL_THROW("A newPtr_ is still active");
      }

      #ifdef UTIL_MPI
      if (cacheSize_ == cacheCapacity_) {
          bool isComplete = false;
          int  source = 0;
          bufferPtr_->endSendBlock(isComplete);
          bufferPtr_->bcast(domainPtr_->communicator(), source);
          nSentTotal_ += cacheSize_;
          bufferPtr_->clearSendBuffer();
          bufferPtr_->beginSendBlock(Buffer::GROUP2 + N - 2);
          cacheSize_ = 0;
      }
      #endif
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
      if (groupStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (domainPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (!domainPtr_->isInitialized()) {
         UTIL_THROW("Domain is not initialized");
      }
      if (domainPtr_->gridRank() != 0) {
         UTIL_THROW("GroupDistributor::add called on slave node");
      }
      if (cache_.capacity() <= 0) {
         UTIL_THROW("GroupDistributor cache is not allocated");
      }
      if (atomStoragePtr_->nGhost() != 0) {
         UTIL_THROW("AtomStorage has ghosts");
      }
      if (newPtr_ == 0) {
         UTIL_THROW("newPtr is null on entry to add()");
      }

      // If group has at least one atoms on master, add to groupStorage.
      int nAtom = atomStoragePtr_->findGroupAtoms(*newPtr_);
      if (nAtom > 0) {
         Group<N>* ptr = groupStoragePtr_->newPtr();
         *ptr = *newPtr_;
         groupStoragePtr_->add();
         nAtomRecv_ += nAtom;
      }
      newPtr_->pack(*bufferPtr_);
 
      // Nullify newPtr_ to release for reuse.
      newPtr_ = 0;
      ++cacheSize_;

   }

   /*
   * Send any groups that have not been sent previously.
   *
   * This method must be called only by the master processor.
   */
   template <int N>
   void GroupDistributor<N>::send()
   {

      // Preconditions
      if (atomStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (groupStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (domainPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (bufferPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("GroupDistributor is not allocated");
      }
      if (domainPtr_->gridRank() != 0) {
         UTIL_THROW("GroupDistributor::send called on slave node");
      }
      if (newPtr_ != 0) {
         UTIL_THROW("A newPtr_ is still active");
      }

      #ifdef UTIL_MPI
      bool isComplete = true;
      int source = 0;
      bufferPtr_->endSendBlock(isComplete);
      bufferPtr_->bcast(domainPtr_->communicator(), source);
      nSentTotal_ += cacheSize_;
      bufferPtr_->clearSendBuffer();

      validate();
      groupStoragePtr_->unsetNTotal();
      groupStoragePtr_->computeNTotal(domainPtr_->communicator());
      if (groupStoragePtr_->nTotal() != nSentTotal_) {
         UTIL_THROW("Number of groups not equal number sent");
      }
      groupStoragePtr_->isValid(*atomStoragePtr_, domainPtr_->communicator(), 
                                false);
      #else
      groupStoragePtr_->unsetNTotal();
      groupStoragePtr_->computeNTotal();
      if (groupStoragePtr_->nTotal() != nSentTotal_) {
         UTIL_THROW("Number of groups not equal number sent");
      }
      groupStoragePtr_->isValid(*atomStoragePtr_, false);
      #endif

      cacheSize_ = 0;
      newPtr_ = 0;
      nSentTotal_ = 0;
      nAtomRecv_ = 0;
   }

   /*
   * Receive all atoms sent by the master processor.
   *
   * Called by all processors except the master.
   */ 
   template <int N>
   void GroupDistributor<N>::receive()
   {
      #ifdef UTIL_MPI
      // Preconditions 
      if (atomStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (groupStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (domainPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (bufferPtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (domainPtr_->gridRank() == 0) {
         UTIL_THROW("GroupDistributor::receive called on master node");
      }
      if (atomStoragePtr_->nGhost() != 0) {
         UTIL_THROW("AtomStorage has ghosts");
      }

      Group<N>* ptr;
      int  nAtom;
      bool isComplete = false;
      const int source = 0;

      nAtomRecv_ = 0;
      while (!isComplete) {

         // Receive broadcast
         bufferPtr_->bcast(domainPtr_->communicator(), source);
      
         // Unpack the buffer, set pointers to atoms
         isComplete = bufferPtr_->beginRecvBlock();
         while (bufferPtr_->recvSize() > 0) {
            ptr = groupStoragePtr_->newPtr();
            ptr->unpack(*bufferPtr_);
            nAtom = atomStoragePtr_->findGroupAtoms(*ptr);
            if (nAtom > 0) {
               groupStoragePtr_->add();
               nAtomRecv_ += nAtom;
            } else if (nAtom == 0) {
               groupStoragePtr_->returnPtr();
            } else {
               UTIL_THROW("Invalid return value from findGroupAtoms");
            }
         }

      }

      // Validate Data
      validate(); // Check number of local atoms in groups
      groupStoragePtr_->unsetNTotal();
      groupStoragePtr_->computeNTotal(domainPtr_->communicator());
      groupStoragePtr_->isValid(*atomStoragePtr_, domainPtr_->communicator(),
                                false);

      nAtomRecv_ = 0;
      #endif
   }

   #ifdef UTIL_MPI
   /**
   * Check number of groups sent and received.
   */
   template <int N>
   void GroupDistributor<N>::validate() 
   {
      int nAtomRecvTot;
      const int source = 0;
      domainPtr_->communicator()
             .Reduce(&nAtomRecv_, &nAtomRecvTot, 1, MPI::INT, MPI::SUM, source);
      if (domainPtr_->gridRank() == 0) {
         if (nAtomRecvTot != nSentTotal_*N) {
             Log::file() << "nSentTotal_*N = " << nSentTotal_*N << std::endl;
             Log::file() << "nAtomRecvTot  = " << nAtomRecvTot  << std::endl;
             UTIL_THROW("Discrepancy in number of local atoms in groups");
         }
      }
   }
   #endif

}
#endif
