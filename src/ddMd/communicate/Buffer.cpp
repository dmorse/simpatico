#ifndef DDMD_BUFFER_CPP
#define DDMD_BUFFER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Buffer.h"
#include "Domain.h"
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Group.h>
#include <util/format/Int.h>
#include <util/mpi/MpiLogger.h>

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
   {}

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

   /**
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
      sendBufferBegin_ = new char[bufferCapacity_];
      if (sendBufferBegin_ == 0) {
         UTIL_THROW("Error: memory not allocated for send buffer");
      }
      sendBufferEnd_ = sendBufferBegin_ + bufferCapacity_;

      // Allocate memory for the receive buffer
      recvBufferBegin_ = new char[bufferCapacity_];
      if (recvBufferBegin_ == 0) {
         UTIL_THROW("Error: memory not allocated for recv buffer");
      }
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

   /*
   * Send and receive buffer.
   */
   void Buffer::sendRecv(MPI::Intracomm& comm, int source, int dest)
   {

      MPI::Request request[2];
      int  sendBytes = 0;
      int  myRank    = comm.Get_rank();
      int  comm_size = comm.Get_size();

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
      request[0] = comm.Irecv(recvBufferBegin_, bufferCapacity_ , 
                              MPI::CHAR, source, 5);

      // Start nonblocking send.
      sendBytes = sendPtr_ - sendBufferBegin_;
      request[1] = comm.Isend(sendBufferBegin_, sendBytes , MPI::CHAR, dest, 5);

      // Wait for completion of receive.
      request[0].Wait();
      recvPtr_ = recvBufferBegin_;

      // Wait for completion of send.
      request[1].Wait();

      // Update statistics.
      if (sendBytes > maxSendLocal_) {
         maxSendLocal_ = sendBytes;
      }
   }

   /*
   * Send a buffer.
   */
   void Buffer::send(MPI::Intracomm& comm, int dest)
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
         UTIL_THROW("Source and destination identical");
      }

      sendBytes = sendPtr_ - sendBufferBegin_;
      request = comm.Isend(sendBufferBegin_, sendBytes, MPI::CHAR, dest, 5);
      request.Wait();

      // Update statistics.
      if (sendBytes > maxSendLocal_) {
         maxSendLocal_ = sendBytes;
      }
   }

   /*
   * Receive a buffer.
   */
   void Buffer::recv(MPI::Intracomm& comm, int source)
   {
      MPI::Request request;
      int  myRank     = comm.Get_rank();
      int  comm_size  = comm.Get_size();

      // Preconditons
      if (source > comm_size - 1 || source < 0) {
         UTIL_THROW("Source rank out of bounds");
      }
      if (source == myRank) {
         UTIL_THROW("Source and destination identical");
      }

      request = comm.Irecv(recvBufferBegin_, bufferCapacity_, 
                           MPI::CHAR, source, 5);
      request.Wait();
      recvType_ = NONE;
      recvPtr_ = recvBufferBegin_;
   }

   /*
   * Broadcast a buffer.
   */
   void Buffer::bcast(MPI::Intracomm& comm, int source)
   {
      int comm_size = comm.Get_size();
      int myRank = comm.Get_rank();
      if (source > comm_size - 1 || source < 0) {
         UTIL_THROW("Source rank out of bounds");
      }

      int sendBytes;
      if (myRank == source) {
         sendBytes = sendPtr_ - sendBufferBegin_;
         comm.Bcast(&sendBytes, 1, MPI::INT, source);
         comm.Bcast(sendBufferBegin_, sendBytes, MPI::CHAR, source);
         sendPtr_ = sendBufferBegin_;
         sendType_ = NONE;
      } else {
         comm.Bcast(&sendBytes, 1, MPI::INT, source);
         comm.Bcast(recvBufferBegin_, sendBytes, MPI::CHAR, source);
         recvPtr_ = recvBufferBegin_;
         recvType_ = NONE;
      }
      if (sendBytes > maxSendLocal_) {
         maxSendLocal_ = sendBytes;
      }

   }

   /*
   * Compute maximum message size among all processors.
   */
   #ifdef UTIL_MPI
   void Buffer::computeStatistics(MPI::Intracomm& comm)
   #else
   void Buffer::computeStatistics()
   #endif
   {
      #ifdef UTIL_MPI
      int globalSendMax;
      comm.Allreduce(&maxSendLocal_, &globalSendMax, 1, MPI::INT, MPI::MAX);
      maxSend_.set(globalSendMax);
      #else
      maxSend_.set(maxSendLocal_);
      #endif
   
   }

   /**
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
   * Number of items currently in data send block. 
   */
   int Buffer::sendSize() const
   {  return sendSize_; }

   /*
   * Number of unread items currently in data receive block. 
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
#endif
