#ifndef PACKED_DATA_CPP
#define PACKED_DATA_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PackedData.h"

namespace Util
{

   /*
   * Constructor.
   */
   PackedData::PackedData()
    : begin_(0),
      endAllocated_(0),
      endPacked_(0),
      cursor_(0),
      capacity_(0),
      status_(NONE)
   {}

   /*
   * Destructor.
   */
   PackedData::~PackedData()
   {  if (begin_) delete begin_; }

   /*
   * Allocate a block of memory.
   */
   void PackedData::allocate(int capacity)
   {
      if (capacity < 0) {
         UTIL_THROW("Negative capacity");
      }
      begin_ = new Byte[capacity];
      endAllocated_ = begin_ + capacity; 
      endPacked_= begin_;
      cursor_ = begin_;
      capacity_ = capacity;
      status_ = NONE;
   }

   /*
   * Prepare to pack data.
   */
   void PackedData::beginPacking()
   {
      cursor_ = begin_; 
      endPacked_ = begin_; 
      status_ = PACKING;
   }

   /*
   * Prepare to unpack data.
   */
   void PackedData::beginUnpacking()
   { 
      endPacked_ = cursor_; 
      cursor_ = begin_; 
      status_ = UNPACKING;
   }

   #if 0
   #ifdef UTIL_MPI
   /*
   * Send a block.
   */
   void PackedData::send(MPI::Intracomm& comm, int dest)
   {
      MPI::Request request;
      int  sendBytes = 0;
      int  comm_size = comm.Get_size();
      int  myRank = comm.Get_rank();

      // Preconditions
      if (dest > comm_size - 1 || dest < 0) {
         UTIL_THROW("Destination rank out of bounds");
      }
      if (dest == myRank) {
         UTIL_THROW("Source and desination identical");
      }

      sendBytes = cursor_ - begin_;
      request = comm.Isend(begin_, sendBytes, MPI::UNSIGNED_CHAR, dest, 5);
      request.Wait();

   }

   /*
   * Receive a block.
   */
   void PackedData::recv(MPI::Intracomm& comm, int source)
   {
      MPI::Request request;
      int  myRank     = comm.Get_rank();
      int  comm_size  = comm.Get_size();

      // Preconditons
      if (source > comm_size - 1 || source < 0) {
         UTIL_THROW("Source rank out of bounds");
      }
      if (source == myRank) {
         UTIL_THROW("Source and desination identical");
      }

      request = comm.Irecv(begin_, capacity_, MPI::UNSIGNED_CHAR, source, 5);
      request.Wait();
      cursor_ = begin_;

   }
   #endif
   #endif

}
#endif
