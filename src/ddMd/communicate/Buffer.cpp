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
    : Util::ParamComposite(),
      #ifdef UTIL_MPI
      sendBufferBegin_(0),
      recvBufferBegin_(0),
      sendBufferEnd_(0),
      recvBufferEnd_(0),
      sendBlockBegin_(0),
      recvBlockBegin_(0),
      sendPtr_(0),
      recvPtr_(0),
      bufferCapacity_(-1),
      dataCapacity_(-1),
      sendSize_(0),
      recvSize_(0),
      #endif
      atomCapacity_(-1),
      ghostCapacity_(-1),
      maxSendLocal_(0),
      isInitialized_(false)
   { setClassName("Buffer"); }

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
      int atomDataSize  = atomCapacity_*atomSize();
      int ghostDataSize = ghostCapacity_*ghostSize();
      if (atomDataSize > ghostDataSize) {
          dataCapacity_ = atomDataSize;
          ghostCapacity_ = dataCapacity_/ghostSize();
      } else {
          dataCapacity_ = ghostDataSize;
          atomCapacity_ = dataCapacity_/atomSize();
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
   void Buffer::beginSendBlock(BlockDataType sendType, int sendGroupSize)
   {
      if (sendSize_ != 0) {
         UTIL_THROW("Error: previous send block not finalized");
      }

      // Set number of atoms currently in send block to zero.
      sendSize_ = 0;

      // Data type to be sent.
      sendType_ = sendType;
      sendGroupSize_ = sendGroupSize; // 0 by default, for atom or ghost.

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
      // Add passport to the beginning of the block:
      // Pack sendSize_, sendType_ and isComplete.
      int* sendBuffPtr = (int *)sendBlockBegin_;
      *sendBuffPtr = sendSize_;
      ++sendBuffPtr;
      *sendBuffPtr = (int) sendType_;
      ++sendBuffPtr;
      *sendBuffPtr = sendGroupSize_;
      ++sendBuffPtr;
      *sendBuffPtr = (int) isComplete;

      // Clear variables associated with the sent block.
      sendBlockBegin_ = 0;
      sendSize_ = 0;
      sendType_ = NONE;
      sendGroupSize_ = 0;

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

      // Extract the number packed items and item type from the receive buffer
      int* recvBuffPtr = (int *)recvPtr_;
      recvSize_  = *recvBuffPtr;
      recvType_  = *(recvBuffPtr + 1);
      recvGroupSize_ = *(recvBuffPtr + 2);
      bool isComplete = (bool) *(recvBuffPtr + 3);
      recvBuffPtr += 4;

      // Set recvPtr to beginning of first item to be unpacked
      recvPtr_ = (char *)recvBuffPtr;

      return isComplete;
   }

   /*
   * Pack a local Atom for exchange of ownership.
   */
   void Buffer::packAtom(Atom& atom)
   {
      if (sendType_ != ATOM) {
         UTIL_THROW("Send type is not ATOM");
      }
      if (sendSize_ >= atomCapacity_) {
         UTIL_THROW("Attempt to overpack send buffer");
      }

      pack<int>(atom.id());
      pack<int>(atom.typeId());
      pack<Vector>(atom.position());
      pack<Vector>(atom.velocity());
      pack<unsigned int>(atom.plan().flags());

      // Pack Mask
      Mask& mask = atom.mask();
      int size = mask.size();
      pack<int>(size);
      for (int j = 0; j < size; ++j) {
         pack<int>(mask[j]);
      }

      //Increment number of atoms in send buffer by 1
      ++sendSize_;
   }

   /*
   * Receive ownership of an Atom.
   */
   void Buffer::unpackAtom(Atom& atom)
   {
      if (recvType_ != (int)ATOM) {
         UTIL_THROW("Receive type is not ATOM");
      }
      if (recvSize_ <= 0) {
         UTIL_THROW("Attempt to unpack empty receive buffer");
      }

      int i;

      unpack(i);
      atom.setId(i);

      unpack(i);
      atom.setTypeId(i);

      unpack<Vector>(atom.position());
      unpack<Vector>(atom.velocity());

      // Unpack communication plan
      unsigned int ui;
      unpack(ui);
      atom.plan().setFlags(ui);

      // Unpack atom Mask
      Mask& mask = atom.mask();
      mask.clear();
      int size;
      unpack(size);
      for (int j = 0; j < size; ++j) {
         unpack(i);
         mask.append(i);
      }
      assert(mask.size() == size);

      // Decrement number of atoms in recv buffer by 1
      recvSize_--;

   }

   /*
   * Return size of packed Atom, in bytes.
   */
   int Buffer::atomSize()
   {  
      int size = 0;
      size += 2*sizeof(int); 
      size += 2*sizeof(Vector); 
      size += sizeof(unsigned int);
      size += sizeof(int);
      size += Mask::Capacity*sizeof(int); 
      return size;
   }

   /*
   * Pack data required for a ghost Atom for sending.
   */
   void Buffer::packGhost(Atom& atom)
   {
      // Preconditions
      if (sendType_ != GHOST) {
         UTIL_THROW("Send type is not GHOST");
      }
      if (sendSize_ >= ghostCapacity_) {
         UTIL_THROW("Attempt to overpack send buffer");
      }

      Vector pos;
      pack<int>(atom.id());
      pack<int>(atom.typeId());
      pack<Vector>(atom.position());
      pack<unsigned int>(atom.plan().flags());

      //Increment number of atoms in send buffer by 1
      sendSize_++;

   }

   /*
   * Unpack data required for a ghost Atom.
   */
   void Buffer::unpackGhost(Atom& atom)
   {
      if (recvType_ != (int)GHOST) {
         UTIL_THROW("Receive type is not GHOST");
      }
      if (recvSize_ <= 0) {
         UTIL_THROW("Attempt to unpack empty receive buffer");
      }

      int i;
      unpack(i);
      atom.setId(i);
      unpack(i);
      atom.setTypeId(i);
      unpack<Vector>(atom.position());

      unsigned int ui;
      unpack(ui);
      atom.plan().setFlags(ui);
      
      //Decrement number of atoms in recv buffer to be unpacked by 1
      recvSize_--;
   }

   /*
   * Return size of one packed Ghost atom, in bytes (static method).
   */
   int Buffer::ghostSize()
   {  
      int size = 0;
      size += 2*sizeof(int); 
      size += sizeof(Vector); 
      size += sizeof(unsigned int);
      return size;
   }

   /*
   * Packed size of one Group<N> object.
   */
   int Buffer::groupSize(int N)
   {  return (2 + N)*sizeof(int) + sizeof(unsigned int); }

   /*
   * Maximum number of Group<N> objects that can fit buffer.
   */
   int Buffer::groupCapacity(int N) const
   {
      if (dataCapacity_ <= 0) {
         UTIL_THROW("Buffer not allocated");
      }
      return (dataCapacity_/groupSize(N)); 
   }

   /*
   * Pack updates ghost atom position.
   */
   void Buffer::packUpdate(Atom& atom)
   {
      // Preconditions
      if (sendType_ != UPDATE) {
         UTIL_THROW("Send type is not UPDATE");
      }
      if (sendSize_ >= ghostCapacity_) {
         UTIL_THROW("Attempt to overpack buffer: sendSize> >= ghostCapacity_");
      }
      pack<Vector>(atom.position());

      //Increment number of atoms in send buffer by 1
      sendSize_++;
   }

   /**
   * Pack updated ghost atom position.
   */
   void Buffer::unpackUpdate(Atom& atom)
   {
      if (recvType_ != (int)UPDATE) {
         UTIL_THROW("Receive type is not UPDATE");
      }
      if (recvSize_ <= 0) {
         UTIL_THROW("Attempt to unpack empty receive buffer");
      }
      unpack<Vector>(atom.position());
      
      //Decrement number of atoms in recv buffer to be unpacked by 1
      recvSize_--;
   }

   /*
   * Pack ghost atom force.
   */
   void Buffer::packForce(Atom& atom)
   {
      // Preconditions
      if (sendType_ != FORCE) {
         UTIL_THROW("Send type is not FORCE");
      }
      if (sendSize_ >= ghostCapacity_) {
         UTIL_THROW("Attempt to overpack buffer: sendSize> >= ghostCapacity_");
      }
      pack<Vector>(atom.force());

      //Increment number of atom forces in send buffer by 1
      sendSize_++;
   }

   /**
   * Unpack data ghost Atom force, and add to atom on this processor.
   */
   void Buffer::unpackForce(Atom& atom)
   {
      if (recvType_ != (int)FORCE) {
         UTIL_THROW("Receive type is not FORCE");
      }
      if (recvSize_ <= 0) {
         UTIL_THROW("Attempt to unpack empty receive buffer");
      }
      Vector f;
      unpack<Vector>(f);
      atom.force() += f;
      
      //Decrement number of atoms in recv buffer to be unpacked by 1
      recvSize_--;
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
