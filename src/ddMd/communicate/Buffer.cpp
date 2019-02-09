/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Buffer.h"
#include "Domain.h"
#include <util/misc/Memory.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Group.h>
#include <util/format/Int.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   Buffer::Buffer()
    : ParamComposite(),
      #ifdef UTIL_MPI
      sendBufferBegin_(0),
      recvBufferBegin_(0),
      sendBufferEnd_(0),
      recvBufferEnd_(0),
      sendBlockBegin_(0),
      recvBlockBegin_(0),
      recvBlockEnd_(0),
      sendPtr_(0),
      recvPtr_(0),
      bufferCapacity_(-1),
      dataCapacity_(-1),
      sendSize_(0),
      recvSize_(0),
      recvType_(NONE),
      #endif
      atomCapacity_(-1),
      ghostCapacity_(-1),
      maxSendLocal_(0),
      isInitialized_(false)
   {  setClassName("Buffer"); }

   /*
   * Destructor.
   */
   Buffer::~Buffer()
   {
      if (sendBufferBegin_) {
         Memory::deallocate<char>(sendBufferBegin_, bufferCapacity_);
      }
      if (recvBufferBegin_) {
         Memory::deallocate<char>(recvBufferBegin_, bufferCapacity_);
      }
   }

   /*
   * Allocate send and recv buffers.
   */
   void Buffer::allocate(int atomCapacity, int ghostCapacity)
   {
      //Preconditions
      if (atomCapacity < 0) {
         UTIL_THROW("Negative atomCapacity");
      }
      if (ghostCapacity < 0) {
         UTIL_THROW("Negative ghostCapacity");
      }

      atomCapacity_  = atomCapacity;
      ghostCapacity_ = ghostCapacity;

      #ifdef UTIL_MPI
      // Do actual allocation
      allocate();
      #endif

      isInitialized_ = true;
   }

   /*
   * Read capacities, and allocate buffers.
   */
   void Buffer::readParameters(std::istream& in)
   {

      // Read parameters
      read<int>(in, "atomCapacity",  atomCapacity_);
      read<int>(in, "ghostCapacity", ghostCapacity_);

      //Preconditions
      if (atomCapacity_ < 0) {
         UTIL_THROW("Negative atomCapacity");
      }
      if (ghostCapacity_ < 0) {
         UTIL_THROW("Negative ghostCapacity");
      }

      #ifdef UTIL_MPI
      // Do actual allocation
      allocate();
      #endif

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void Buffer::loadParameters(Serializable::IArchive &ar)
   {
      // Read parameters
      loadParameter<int>(ar, "atomCapacity",  atomCapacity_);
      loadParameter<int>(ar, "ghostCapacity", ghostCapacity_);

      // Validate data
      if (atomCapacity_ < 0) {
         UTIL_THROW("Negative atomCapacity");
      }
      if (ghostCapacity_ < 0) {
         UTIL_THROW("Negative ghostCapacity");
      }

      #ifdef UTIL_MPI
      // Do actual allocation
      allocate();
      #endif

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void Buffer::save(Serializable::OArchive &ar)
   {
      ar << atomCapacity_;
      ar << ghostCapacity_;
   }

   /*
   * Maximum number of atoms for which space is available.
   */
   int Buffer::atomCapacity() const
   {  return atomCapacity_; }

   /*
   * Maximum number of ghost atoms for which space is available.
   */
   int Buffer::ghostCapacity() const
   {  return ghostCapacity_; }

   /*
   * Has this buffer been initialized?
   */
   bool Buffer::isInitialized() const
   {  return isInitialized_; }

   #ifdef UTIL_MPI
   /*
   * Allocate send and recv buffers (private method).
   *
   * This method uses values of atomCapacity_ and ghostCapacity_ that
   * must have been set previously. It is called by allocate(int, int)
   * and readParameters() to do the actual allocation.
   */
   void Buffer::allocate()
   {

      // Preconditions
      if (atomCapacity_ <= 0) {
         UTIL_THROW("atomCapacity_ must be positive");
      }
      if (ghostCapacity_ <= 0) {
         UTIL_THROW("ghostCapacity_ must be positive");
      }
      if (isAllocated()) {
         UTIL_THROW("Buffer cannot be re-allocated");
      }

      // Capacity in bytes send and receive buffers. This is maximum of
      // the buffer space required by the local atoms and ghost atoms.
      int atomDataSize  = atomCapacity_*Atom::packedAtomSize();
      int ghostDataSize = ghostCapacity_*Atom::packedGhostSize();
      if (atomDataSize > ghostDataSize) {
          dataCapacity_ = atomDataSize;
          ghostCapacity_ = dataCapacity_/Atom::packedGhostSize();
      } else {
          dataCapacity_ = ghostDataSize;
          atomCapacity_ = dataCapacity_/Atom::packedAtomSize();
      }

      // Leave space for a 4 byte header
      bufferCapacity_ += dataCapacity_ + 4 * sizeof(int);

      // Allocate memory for the send buffer
      Memory::allocate<char>(sendBufferBegin_, bufferCapacity_);
      sendBufferEnd_ = sendBufferBegin_ + bufferCapacity_;

      // Allocate memory for the receive buffer
      Memory::allocate<char>(recvBufferBegin_, bufferCapacity_);
      recvBufferEnd_ = recvBufferBegin_ + bufferCapacity_;

      recvPtr_ = recvBufferBegin_;
   }

   /*
   * Clear the send buffer prior to packing, and set the sendType.
   */
   void Buffer::clearSendBuffer()
   {
      sendPtr_ = sendBufferBegin_;
      sendSize_ = 0;
      sendType_ = NONE;
      sendBlockBegin_ = 0;
   }

   /*
   * Clear the send buffer prior to packing, and set the sendType.
   */
   void Buffer::beginSendBlock(int sendType)
   {
      if (sendSize_ != 0) {
         UTIL_THROW("Error: previous send block not finalized");
      }

      // Set number of atoms currently in send block to zero.
      sendSize_ = 0;

      // Data type to be sent.
      sendType_ = sendType;

      // Mark beginning of block.
      sendBlockBegin_ = sendPtr_;

      // Increment sendPtr_ to leave space for 4 integers.
      int* sendBuffPtr = (int *)sendPtr_;
      sendBuffPtr += 4;
      sendPtr_ = (char *)sendBuffPtr;
   }

   /*
   * Finalize data block in buffer. Pack prefix data.
   */
   void Buffer::endSendBlock(bool isComplete)
   {
      int sendBytes = (int)(sendPtr_ - sendBlockBegin_);

      // Add passport to the beginning of the block:
      // Pack sendSize_, sendBytes, sendType_ and isComplete.
      int* sendBuffPtr = (int *)sendBlockBegin_;
      *sendBuffPtr = sendSize_;
      ++sendBuffPtr;
      *sendBuffPtr = sendBytes;
      ++sendBuffPtr;
      *sendBuffPtr = sendType_;
      ++sendBuffPtr;
      *sendBuffPtr = (int) isComplete;

      // Clear variables associated with the sent block.
      sendBlockBegin_ = 0;
      sendSize_ = 0;
      sendType_ = NONE;
   }

   /*
   * Begin receiving block. Extract prefix data.
   */
   bool Buffer::beginRecvBlock()
   {
      // Precondition
      if (recvSize_ != 0) {
         UTIL_THROW("Error: Previous receive block not completely unpacked");
      }

      // Store address of begining of block
      recvBlockBegin_ = recvPtr_;

      // Extract passport data
      int* recvBuffPtr = (int *)recvPtr_;
      recvSize_  = *recvBuffPtr;
      int recvBytes = *(recvBuffPtr + 1);
      recvType_  = *(recvBuffPtr + 2);
      bool isComplete = (bool) *(recvBuffPtr + 3);

      // Calculate address of expected end of recv block
      recvBlockEnd_ = recvBlockBegin_ + recvBytes;

      // Set recvPtr to beginning of first item to be unpacked
      recvBuffPtr += 4;
      recvPtr_ = (char *)recvBuffPtr;

      return isComplete;
   }

   /*
   * Finalize receiving block, check consistency.
   */
   void Buffer::endRecvBlock()
   {
      if (recvSize_ != 0) {
         UTIL_THROW("Error: Recv counter != 0 at end of block");
      }
      if (recvPtr_ != recvBlockEnd_) {
         UTIL_THROW("Error: Inconsistent recv cursor at end of block");
      }
      recvBlockBegin_ = 0;
      recvBlockEnd_ = 0;
      recvSize_ = 0;
      recvType_ = NONE;
   }

   #ifdef UTIL_MPI
   /*
   * Send and receive buffer.
   */
   void Buffer::sendRecv(MPI_Comm comm, int source, int dest)
   {

      MPI_Request request[2];
      MPI_Status status[2];
      int  sendBytes = 0;
      int  myRank;
      MPI_Comm_rank(comm, &myRank);
      int  comm_size;
      MPI_Comm_size(comm, &comm_size);

      // Preconditions
      if (dest > comm_size - 1 || dest < 0) {
         UTIL_THROW("Destination rank out of bounds");
      }
      if (source > comm_size - 1 || source < 0) {
         UTIL_THROW("Source rank out of bounds");
      }
      if (dest == myRank) {
         UTIL_THROW("Destination and my rank are identical");
      }
      if (source == myRank) {
         UTIL_THROW("Source and my rank are identical");
      }

      // Start nonblocking receive.
      MPI_Irecv(recvBufferBegin_, bufferCapacity_ , MPI_CHAR, source, 5, 
                comm, &request[0]);

      // Start nonblocking send.
      sendBytes = sendPtr_ - sendBufferBegin_;
      MPI_Isend(sendBufferBegin_, sendBytes , MPI_CHAR, dest, 5, 
                comm, &request[1]);

      // Wait for completion of receive.
      //request[0].Wait();
      MPI_Wait(&request[0], &status[0]);
      recvPtr_ = recvBufferBegin_;

      // Wait for completion of send.
      // request[1].Wait();
      MPI_Wait(&request[1], &status[1]);

      // Update statistics.
      if (sendBytes > maxSendLocal_) {
         maxSendLocal_ = sendBytes;
      }
   }

   /*
   * Send a buffer.
   */
   void Buffer::send(MPI_Comm comm, int dest)
   {
      int  myRank;
      MPI_Comm_rank(comm, &myRank);
      int  comm_size;
      MPI_Comm_size(comm, &comm_size);

      // Preconditions
      if (dest > comm_size - 1 || dest < 0) {
         UTIL_THROW("Destination rank out of bounds");
      }
      if (dest == myRank) {
         UTIL_THROW("Source and destination identical");
      }

      int  sendBytes = 0;
      sendBytes = sendPtr_ - sendBufferBegin_;
      MPI_Request request;
      MPI_Isend(sendBufferBegin_, sendBytes, MPI_CHAR, dest, 5, 
                comm, &request);
      MPI_Status status;
      MPI_Wait(&request, &status);

      // Update statistics.
      if (sendBytes > maxSendLocal_) {
         maxSendLocal_ = sendBytes;
      }
   }

   /*
   * Receive a buffer.
   */
   void Buffer::recv(MPI_Comm comm, int source)
   {
      int myRank;
      MPI_Comm_rank(comm, &myRank);
      int comm_size;
      MPI_Comm_size(comm, &comm_size);

      // Preconditons
      if (source > comm_size - 1 || source < 0) {
         UTIL_THROW("Source rank out of bounds");
      }
      if (source == myRank) {
         UTIL_THROW("Source and destination identical");
      }

      MPI_Request request;
      MPI_Irecv(recvBufferBegin_, bufferCapacity_, MPI_CHAR, source, 5, 
                comm, &request);
      MPI_Status status;
      MPI_Wait(&request, &status);
      recvType_ = NONE;
      recvPtr_ = recvBufferBegin_;
   }

   /*
   * Broadcast a buffer.
   */
   void Buffer::bcast(MPI_Comm comm, int source)
   {
      int myRank;
      MPI_Comm_rank(comm, &myRank);
      int comm_size;
      MPI_Comm_size(comm, &comm_size);
      if (source > comm_size - 1 || source < 0) {
         UTIL_THROW("Source rank out of bounds");
      }

      int sendBytes;
      if (myRank == source) {
         sendBytes = sendPtr_ - sendBufferBegin_;
         MPI_Bcast(&sendBytes, 1, MPI_INT, source, comm);
         MPI_Bcast(sendBufferBegin_, sendBytes, MPI_CHAR, source, comm);
         sendPtr_ = sendBufferBegin_;
         sendType_ = NONE;
      } else {
         MPI_Bcast(&sendBytes, 1, MPI_INT, source, comm);
         MPI_Bcast(recvBufferBegin_, sendBytes, MPI_CHAR, source, comm);
         recvPtr_ = recvBufferBegin_;
         recvType_ = NONE;
      }
      if (sendBytes > maxSendLocal_) {
         maxSendLocal_ = sendBytes;
      }

   }
   #endif

   /*
   * Compute maximum message size among all processors.
   */
   #ifdef UTIL_MPI
   void Buffer::computeStatistics(MPI_Comm comm)
   #else
   void Buffer::computeStatistics()
   #endif
   {
      #ifdef UTIL_MPI
      int globalSendMax;
      MPI_Allreduce(&maxSendLocal_, &globalSendMax, 1, MPI_INT, MPI_MAX, comm);
      maxSend_.set(globalSendMax);
      #else
      maxSend_.set(maxSendLocal_);
      #endif
   }

   /*
   * Clear any accumulated usage statistics.
   */
   void Buffer::clearStatistics()
   {
      maxSendLocal_ = 0;
      maxSend_.unset();
   }

   /*
   * Output statistics.
   */
   void Buffer::outputStatistics(std::ostream& out)
   {

      out << std::endl;
      out << "Buffer" << std::endl;
      out << "sendBytes: max, capacity "
          << Int(maxSend_.value(), 10)
          << Int(bufferCapacity_, 10)
          << std::endl;
   }

   /*
   * Number of items packed thus far in current data send block.
   */
   int Buffer::sendSize() const
   {  return sendSize_; }

   /*
   * Number of unread items left in the current receive block.
   */
   int Buffer::recvSize() const
   {  return recvSize_; }

   /*
   * Has this buffer been allocated?
   */
   bool Buffer::isAllocated() const
   {  return (bufferCapacity_ > 0); }
   #endif

}
