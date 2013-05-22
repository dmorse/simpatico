#ifndef DDMD_BUFFER_H
#define DDMD_BUFFER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <util/misc/Setable.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   class Atom;
   template <int N> class Group;

   /**
   * Buffer for sending atoms and groups between processors.
   *
   * \ingroup DdMd_Communicate_Module
   */
   class Buffer: public ParamComposite 
   {

   public:

      enum BlockDataType {NONE, ATOM, GHOST, UPDATE, FORCE, GROUP};

      /**
      * Constructor.
      */
      Buffer();

      virtual ~Buffer();

      /**
      * Read capacities and allocate buffers.
      *
      * Read parameters atomCapacity and ghostCapacity (see allocate)
      * and allocate buffers.
      */
      void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /**
      * Allocate send and recv buffers.
      *
      * Calculate amounts of memory required to accommodate specified
      * numbers of local atoms and ghosts, and allocate buffers using
      * the larger value.
      *
      * \param atomCapacity max expected number of local atoms.
      * \param ghostCapacity max expected number of ghost atoms.
      */
      void allocate(int atomCapacity, int ghostCapacity);

      /// \name Send Buffer Management
      //@{
      
      /**
      * Clear the send buffer and sendType.
      */
      void clearSendBuffer();

      /**
      * Clear the send buffer prior to packing and set atomtype.
      *
      * Sets sendSize() to zero, the send pointer to the beginning of
      * the send buffer, and sendType to the type of data to be sent.
      * If sendType = Group, also set sendGroupSize.
      */
      void beginSendBlock(BlockDataType sendType, int sendGroupSize = 0);

      /**
      * Function template for packing one variable into the send buffer.
      *
      * This is designed for primitive C variables. It will work on any
      * plain old data type T for which the assignment (=) operator does
      * a straight bitwise copy from the variable to the buffer.
      *
      * \param data variable to be packed
      */
      template <typename T>
      void pack(const T& data);

      /**
      * Increment sendSize counter after packing an item.
      */
      void incrementSendSize();

      /**
      * Finalize a block in the send buffer.
      *
      * \param isComplete false if data block is incomplete, true otherwise
      */
      void endSendBlock(bool isComplete = true);

      //@}
      /// \name Receive Buffer Management
      //@{
      
      /**
      * Begin to receive a block from the recv buffer, read envelope.
      *
      * \return false if an incomplete block was received, true otherwise
      */
      bool beginRecvBlock();

      /**
      * Function template unpacking one variable from the receive buffer.
      *
      * This is designed for primitive C variables. It will work on any
      * plain old data type T for which the assignment (=) operator does
      * a straight bitwise copy from the buffer to the variable.
      *
      * \param data variable to be unpacked
      */
      template <typename T>
      void unpack(T& data);

      /**
      * Decrement recvSize counter after unpacking an item.
      */
      void decrementRecvSize();

      //@}
      /// \name Interprocessor Communication
      //@{
      
      /**
      * Receive from processor send and send to processor recv.
      *
      * After transmission, clears the send buffer for reuse, by calling 
      * clearSend(), and sets the receive buffer for unpacking, by setting 
      * recvSize() to the total number of atoms recieved and setting the 
      * receive pointer to first atom to be unpacked. 
      *
      * Throws an Exception if atomType() == NONE or sendSize() == 0.
      *
      * \param comm   MPI communicator object
      * \param source MPI rank of processor from which data is sent
      * \param dest   MPI rank of processor to which data is sent
      */
      void sendRecv(MPI::Intracomm& comm, int source, int dest);

      /**
      * Send a complete buffer.
      *
      * Upon sending buffer, clears the send buffer for reuse by calling
      * clearSend().
      *
      * \param comm  MPI communicator
      * \param dest  rank of processor to which data is sent
      */
      void send(MPI::Intracomm& comm, int dest);

      /**
      * Receive a buffer.
      *
      * \param comm   MPI communicator
      * \param source rank of processor from which data is sent.
      */
      void recv(MPI::Intracomm& comm, int source);

      /**
      * Broadcast a buffer.
      *
      * \param comm   MPI communicator.
      * \param source rank of source processor.
      */
      void bcast(MPI::Intracomm& comm, int source);

      //@}
      /// \name Statistics
      //@{
    
      /**
      * Compute statistics (reduce from all processors).
      * 
      * Call on all processors.
      *
      * \param comm MPI communicator
      */
      #ifdef UTIL_MPI
      virtual void computeStatistics(MPI::Intracomm& comm);
      #else
      virtual void computeStatistics();
      #endif

      /**
      * Output statistics.
      *
      * Call on master, after calling computeStatistics on all procs.
      *
      * \param out   output stream
      */
      void outputStatistics(std::ostream& out);

      /**
      * Clear any accumulated usage statistics.
      */
      void clearStatistics();
      
      //@}
      /// \name Accessors
      //@{
      
      /**
      * Number of items in current send block.
      */
      int sendSize() const;

      /**
      * Number of unread items in current recv block.
      */
      int recvSize() const;

      /**
      * Has memory been allocated for this Buffer?
      */
      bool isAllocated() const;

      /**
      * Has this Buffer been initialized?
      */
      bool isInitialized() const;

      /**
      * Maximum number of atoms for which space is available.
      */
      int atomCapacity() const;

      /**
      * Maximum number of ghost atoms for which space is available.
      */
      int ghostCapacity() const;

      /**
      * Maximum number of group<N> objects for which space is available.
      */
      template <int N>
      int groupCapacity() const;

      //@}

   private:

      /// Pointer to send buffer.
      char* sendBufferBegin_;

      /// Pointer to recv buffer.
      char* recvBufferBegin_;

      /// End of allocated send Buffer (one char past end).
      char* sendBufferEnd_;

      /// End of allocated send Buffer (one char past end).
      char* recvBufferEnd_;

      /// Address one past end of the packed portion of the send buffer.
      char* sendBlockBegin_;

      /// Address one past end of the unpacked portion of the send buffer.
      char* recvBlockBegin_;

      /// Address one past end of the packed portion of the send buffer.
      char* sendPtr_;

      /// Address one past end of the unpacked portion of the send buffer.
      char* recvPtr_;

      /// Allocated capacity of send and recv buffers, in bytes.
      int bufferCapacity_;

      /// Capacity of buffers, in bytes, without 4 int envelope.
      int dataCapacity_;

      /// Number of atoms or ghosts currently in send buffer.
      int sendSize_;

      /// Number of unread atoms or ghosts currently in receive buffer.
      int recvSize_;

      /// Number of atoms in group (or 0 if not a Group).
      int sendGroupSize_;

      /// Type of atom being sent = NONE, ATOM, GHOST, GROUP
      BlockDataType sendType_;

      /// Type of atom being received (BlockDataType cast to int)
      int recvType_;

      /// Number of atoms in a Group type (or 0 if not a Group).
      int recvGroupSize_;

      /// Maximum number of local atoms in buffer.
      int atomCapacity_;

      /// Maximum number of ghost atoms in buffer.
      int ghostCapacity_;

      /// Maximum size used for send buffer on this processor, in bytes.
      int maxSendLocal_;

      /// Has this buffer been initialized ?
      bool isInitialized_;

      /// Maximum size used for send buffers on any processor, in bytes.
      Setable<int> maxSend_;

      /*
      * Allocate send and recv buffers, using preset capacities.
      */
      void allocate();

   };

   /*
   * Pack an object of type T into send buffer.
   */
   template <typename T>
   inline void Buffer::pack(const T& data)
   {
      if (sendPtr_ + sizeof(data) > sendBufferEnd_) {
         UTIL_THROW("Attempted write past end of send buffer");
      }
      T* ptr = (T *)sendPtr_;
      *ptr = data;
      ++ptr;
      sendPtr_ = (char *)ptr;
   }

   /*
   * Unpack an object of type T from recvBuffer.
   */
   template <typename T>
   inline void Buffer::unpack(T& data)
   {
      if (recvPtr_ + sizeof(data) > recvBufferEnd_) {
         UTIL_THROW("Attempted read past end of recv buffer");
      }
      T* ptr = (T *)recvPtr_;
      data = *ptr;
      ++ptr;
      recvPtr_ = (char *)ptr;
   }

   /*
   * Maximum number of Group<N> objects that can fit buffer.
   */
   template <int N>
   int Buffer::groupCapacity() const
   {
      if (dataCapacity_ <= 0) {
         UTIL_THROW("Buffer not allocated");
      }
      return (dataCapacity_/Group<N>::packedSize()); 
   }

   /*
   * Increment sendSize counter after packing an item.
   */
   inline void Buffer::incrementSendSize()
   { ++sendSize_; }

   /*
   * Decrement sendSize counter after receiving an item.
   */
   inline void Buffer::decrementRecvSize()
   { --recvSize_; }

}
#endif
