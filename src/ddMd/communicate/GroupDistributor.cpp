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
    : cache_()
      newPtr_(0),
      atomStoragePtr_(0),
      groupStoragePtr_(0),
      cacheCapacity_(0),
      cacheSize_(0)
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
                                       AtomStorage& atomStorage)
   {
      groupStoragePtr_ = &groupStorage;
      atomStoragePtr_  = &atomStorage;
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

      if (cacheSize_ == cacheCapacity_) {
        
          #if UTIL_MPI
          // Broadcast cache
          #endif
          cacheSize_ = 0;

      }
      return &cache_[i];
   }

   /*
   * Add an atom to the list to be sent.
   */ 
   template <int N>
   void GroupDistributor<N>::add();
   {
      // Preconditions
      if (atomStoragePtr_ == 0) {
         UTIL_THROW("GroupDistributor is not initialized");
      }
      if (cacheCapacity_ <= 0) {
         UTIL_THROW("GroupDistributor is not allocated");
      }

      // Decide if this group is owned by the amster
      // If not owned by the master, queue this atom for sending.
      if (ownedByMaster) {

         Group<N>* ptr  = storage.newPtr();
         *ptr = *newPtr_;
         storage.addNewGroup();

      #ifndef UTIL_MPI
      } else {
         UTIL_THROW("Group not owned by master but UTIL_MPI not defined");
      #endif
      }
    
      // Nullify newPtr_ to release for reuse.
      newPtr_ = 0;

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
         // Broadcast cache;
      }

      // Postconditions
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
   void GroupDistributor<N>::receive()
   {

      #ifdef UTIL_MPI

      // Preconditions

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
            // Do I own any of the atoms in this group ? 
            storage.add();
         }

      }
      #endif

   }

}
#endif
