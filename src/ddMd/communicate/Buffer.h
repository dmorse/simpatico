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

   template <int N> class Group;

   /**
   * Buffer for sending blocks of data between processors.
   *
   * A Buffer manages two blocks of memory blocks, a send buffer
   * and a receive buffer, that are used to communicate data
   * between processors. The class provides a simple interace
   * for: (1) Packing data into the send buffer on a processor
   * from which it will be sent, (2) Unpacking data from the 
   * receive buffer on the processor that receives a buffer, 
   * and (3) sending, receiving, and broadcasting the data in
   * these two buffers. The functions that send and receive
   * data are relatively simple wrappers around MPI functions.
   *
   * A send or receive buffer may contain one or more blocks
   * of data. Each block contains a sequence of items of the 
   * same type. The expected types of data are enumerated by
   * the public enum Buffer::BlockDataType. Each item in a 
   * block may contain an Atom packed for exchange of ownerhsip
   * (ATOM) of ownership, the position, id and type of a ghost 
   * (GHOST) atom, an update of a ghost position (UPDATE), a
   * force vector for use in reverse update (FORCE), or any
   * of several types of covalent Group (GROUP2, GROUP3, and
   * GROUP4). 
   *
   * The DdMd::Atom  class and DdMd::Group class template
   * provide functions to pack and unpack individual items
   * (e.g., atoms or groups). The pack and unpack functions 
   * in these classes are implemented using the primitive 
   * pack() and unpack() function templates, which allow the
   * user to pack and unpack a single primitive C variable 
   * into a heterogeneous buffer. 
   *
   * \ingroup DdMd_Communicate_Module
   */
   class Buffer: public ParamComposite 
   {

   public:

      /**
      * Enumeration of types of data to be sent in blocks. 
      */
      enum BlockDataType 
           {NONE, ATOM, GHOST, UPDATE, FORCE, GROUP2, GROUP3, GROUP4};

      /**
      * Constructor.
      */
      Buffer();

      /**
      * Destructor.
      */
      virtual ~Buffer();

      /**
      * Read capacity and allocate buffers.
      *
      * \param in input parameter stream
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

      /** \name Send Buffer Management
      *
      * To pack a block of data into the send buffer.
      *
      * - Call beginSendBlock(int) at the beginning of the block.
      *
      * - Call the appropriate method of DdMd::Atom or DdMd::Group
      *   once for each item to be packed, within a loop over items 
      *   to be sent.
      *
      * - After packing each item, call incrementSendSize() to
      *   increment a counter of the number of items in the block.
      *
      * - Call endSendBlock(bool) at the end of the block, to
      *   finalize the block. 
      *
      * A complete send buffer may contain one or more such blocks
      * of data, e.g., it may contain a block of atoms followed by
      * several blocks of groups.
      */
      //@{
      
      /**
      * Clear the entire send buffer.
      */
      void clearSendBuffer();

      /**
      * Initialize the block atomtype.
      *
      * Sets sendSize() to zero and sets the sendType.
      *
      * \param sendType BlockDataType value for type of data to be sent
      */
      void beginSendBlock(int sendType);

      /**
      * Function template for packing one variable into the send buffer.
      *
      * This method is used to implement the pack functions provided
      * by the DdMd::Atom and DdMd::Group classes. It is designed to 
      * pack a single primitive C variable into the send buffer. It will 
      * work on any plain old data type T for which the assignment (=) 
      * operator does a straight bitwise copy.
      *
      * \param data variable to be packed
      */
      template <typename T>
      void pack(const T& data);

      /**
      * Increment sendSize counter after packing an item (an Atom or Group).
      */
      void incrementSendSize();

      /**
      * Finalize a block in the send buffer.
      *
      * This method writes a "descriptor" section at the beginning of the block,
      * that describes the associated data block. The descriptor specifies the 
      * length of the block, the type of data, and whether the block is 
      * "complete". 
      *
      * A block should be marked as incomplete iff all of the required data 
      * of the relevant sendtype did not fit into the buffer. This tells the
      * receiving processor to expect one or more other buffers containing 
      * the remaining data of that type.
      *
      * \param isComplete false if data block is incomplete, true otherwise
      */
      void endSendBlock(bool isComplete = true);

      //@}
      /// \name Receive Buffer Management
      //@{
      
      /**
      * Begin to receive a block from the recv buffer.
      *
      * This method reads the descriptor section at the beginning of a
      * block, and sets the receive cursor to be ready to read the first
      * item.
      *
      * \return false if this block is complete, false otherwise.
      */
      bool beginRecvBlock();

      /**
      * Function template unpacking one variable from the receive buffer.
      *
      * This method is used to implement the unpack functions provided
      * by the DdMd::Atom and DdMd::Group classes. It is designed to 
      * unpack a single primitive C variable from the buffer. It will 
      * work on any plain old data type T for which the assignment (=) 
      * operator does a straight bitwise copy.
      *
      * \param data variable into which data should be copied from buffer.
      */
      template <typename T>
      void unpack(T& data);

      /**
      * Decrement recvSize counter after unpacking an item.
      */
      void decrementRecvSize();

      /**
      * Finish processing a block in the recv buffer.
      *
      * This method checks that the recvSize() is zero at the end
      * of a block, and that the receive buffer cursor is at the
      * expected positon, and throws an exception if any suprises
      * occur. It also nullifies the recvSize counter.
      */
      void endRecvBlock();

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

      /// Pointer to beginning of current block in send buffer.
      char* sendBlockBegin_;

      /// Pointer to beginning of current block in recv buffer.
      char* recvBlockBegin_;

      /// Pointer to end of current block in recv buffer.
      char* recvBlockEnd_;

      /// Address one past end of the packed portion of the send buffer.
      char* sendPtr_;

      /// Address one past end of the unpacked portion of the send buffer.
      char* recvPtr_;

      /// Allocated capacity of send and recv buffers, in bytes.
      int bufferCapacity_;

      /// Capacity of buffers, in bytes, without 4 int envelope.
      int dataCapacity_;

      /// Number of items packed thus far into current block of send buffer.
      int sendSize_;

      /// Number of unread items remaining in this block of receive buffer.
      int recvSize_;

      /// Type of item being sent in this block = NONE, ATOM, GHOST, GROUP
      int sendType_;

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

      /// Maximum size used for send buffers on any processor, in bytes.
      Setable<int> maxSend_;

      /// Has this buffer been initialized ?
      bool isInitialized_;

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
