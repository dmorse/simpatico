#ifndef DDMD_BUFFER_H
#define DDMD_BUFFER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <util/misc/Setable.h>          // member
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   template <int N> class Group;

   /**
   * Buffer for interprocessor communication.
   *
   * A Buffer manages two blocks of raw memory, a send buffer
   * and a receive buffer, that are used to communicate data
   * between processors. The class provides a simple interace
   * for: 
   * 
   *   - Packing data into the send buffer on the processor
   *     that sends data (source processor)
   * 
   *   - Unpacking data from the receive buffer on a processor 
   *     that receives the data (destination processor)
   * 
   *   - Sending, receiving, and broadcasting the data in these 
   *     two buffers. 
   *
   * The functions that send, receive and broadcast data are 
   * relatively simple wrappers around MPI functions.
   *
   * The send and receive buffers contain one or more blocks of
   * data. Each block contains a prefix and data segment. The
   * data segment contains a sequence of items of the same data
   * type. The prefix contains information about the data segment,
   * including the length of the buffer, the number of items, 
   * and the type of data.
   *
   * \section Buffer_type_sec Data types
   *
   * Different possible types of data are enumerated by the public 
   * enum Buffer::BlockDataType. The allowed values of this enum are:
   *
   *   - NONE    : Default setting (between blocks)
   *   - ATOM    : a full atom packed for exchange of ownership
   *   - GHOST   : a ghost atom (position, id and typeId)
   *   - UPDATE  : an update of a ghost atom position
   *   - FORCE   : a force vector for use in a reverse update
   *   - GROUP2  : a Group<2> (bond)
   *   - GROUP3  : a Group<3> (angle)
   *   - GROUP4  : a Group<4> (dihedral)
   *   - SPECIAL : specialized nonstandard data type
   *   
   * The standard ATOM, GHOST, UPDATE, FORCE, and GROUP(2,3,4) types
   * are used in the implementation of the Exchanger class.
   *
   * The DdMd::Atom  class and DdMd::Group class template provide
   * provide member functions to pack and unpack individual items 
   * of standard types (atoms and groups). The pack and unpack 
   * functions in these classes each take a Buffer by reference 
   * as an argument. They are implemented using the primitive 
   * Buffer::pack<T>() and Buffer::unpack<T>() member function 
   * templates of the Buffer class, which a user to pack and 
   * unpack a single C variable of type T to or from a Buffer.
   * 
   * The SPECIAL BlockDataType value is a generic label for any
   * specialized, non-standard data type. Code that uses an 
   * Buffer to communicate a specialized data type should provide
   * functions to pack and unpack individual items of the required
   * type. Like the pack and unpack functions for atoms and groups,
   * these functions should be implemented using the Buffer::pack() 
   * and Buffer::unpack() function templates.
   *
   * \section Buffer_pack_sec Packing a block of data
   *
   * Example - packing a block of Group<2> objects:
   * \code
   *
   *    Buffer buffer;
   *    DArray< Group<2> > bonds;
   *
   *    // Clear buffer before packing
   *    buffer.clear();   
   * 
   *    // Pack data block
   *    buffer.beginSendBlock(Buffer::GROUP2);
   *    for (i = 0; i < bonds.size(); ++i) {
   *       bonds[i].pack(buffer);
   *       buffer.incrementSendSize();
   *    }
   *    bool isComplete = true;
   *    endSendBlock(true);
   *
   *    // Pack subsequent blocks, if any.
   *
   *    buffer.send(communicator, destination);
   *
   * \endcode
   * Note:
   *
   *  - Call beginSendBlock() and endSendBlock() before and after the 
   *    loop over data items, respectively.
   *
   *  - Call incrementSendSize() after packing each item to increment
   *    an internal counter of number of items in block. 
   * 
   *  -  The bool isComplete argument of endSendBlock(bool) indicates 
   *     whether the block of data is complete, i.e., whether all data 
   *     of this type is contained in this block (true), or whether
   *     the receiving processor should expect one or more further 
   *     messages containing the rest of this block. This flag is sent
   *     as part of the block prefix.
   *
   * \section Buffer_unpack_sec Unpacking a block of data
   *
   * Example - unpacking a block of Group<2> objects:
   * \code
   *
   *    Buffer buffer;
   *    DArray< Group<2> > bonds;
   *    MPI::Intracomm communicator;
   *    int source;
   *
   *    // Receive the buffer 
   *    buffer.recv(communicator, source);
   *
   *    // Unpack a data block
   *    bool isComplete;
   *    int i = 0;
   *    isComplete = buffer.beginRecvBlock();
   *    while (buffer.recvSize()) {
   *       bonds[i].unpack(buffer);
   *       buffer.decrementRecvSize();
   *       ++i;
   *    }
   *    endRecvBlock();
   *
   *    // Subsequent blocks, if any
   *
   * \endcode
   * Note:
   *
   *  - Call beginRecvBlock() and endRecvBlock() before and after 
   *    the loop over items in the block, respectively. The return 
   *    value of Buffer::recvSize() is equal to the total number 
   *    of data items in the block after beginRecvBlock() returns, 
   *    and must have reached zero when endRecvBlock() is called. 
   *
   *  - The return value of beginRecvBlock() is the value of the
   *    isComplete flag that was passed as an argument to the 
   *    endSendBlock(bool) function on the source processor. 
   *
   *  - Call decrementRecvSize() after unpackingpacking each item 
   *    to decrement the value recvSize counter by one. 
   * 
   *  - This example does not include code to check the isComplete 
   *    flag. More complete code would check this, and either throw
   *    an exception or prepare to receive another message if the
   *    flag was set false.
   *
   * \section Buffer_isComplete_sec Sending a block in several messages
   *
   * The above examples shows a simple case of packing and unpacking
   * a complete block without any explicit checking of whether the 
   * message will fit in the buffer on the source processor or checking
   * of the isComplete flag on the destination. Code that needs to 
   * send blocks that may exceed the buffer size can use the isComplete
   * flag to break the a block of data into in several messages, by 
   * setting isComplete false in all but the last message of the 
   * sequence. 
   * 
   * Here is an example of the code to receive an array of data that 
   * may have been split into several messages:
   * \code
   *
   *    Buffer buffer;
   *    DArray< Group<2> > bonds;
   *    MPI::Intracomm communicator;
   *    int source;
   *
   *    int i = 0;
   *    bool isComplete = false;
   *    while (!isComplete) {
   *       buffer.recv(communicator, source);
   *       isComplete = buffer.beginRecvBlock();
   *       while (buffer.recvSize()) {
   *          bonds[i].unpack(buffer);
   *          buffer.decrementRecvSize();
   *          ++i;
   *       }
   *       endRecvBlock();
   *    }
   *
   * \endcode
   * The code required to avoid overfilling the send buffer on the
   * source processor (not shown), requires knowledge of the packed
   * size of one object. This is provided for Atom and Group objects
   * by the Atom::packedAtomSize(), Atom::packedGhostSize() and
   * Group::packedSize() static functions. 
   *
   * This idiom for dividing large blocks is used in the implementation
   * of all the DdMd Distributor and Collector classes (AtomDistributor, 
   * AtomCollector, GroupDistributor and GroupCollector).
   *
   * \ingroup DdMd_Communicate_Module
   */
   class Buffer: public ParamComposite 
   {

   public:

      /**
      * Enumeration of types of data to be sent in blocks. 
      */
      enum BlockDataType {NONE, ATOM, GHOST, UPDATE, FORCE, 
                          GROUP2, GROUP3, GROUP4, SPECIAL};

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
      * number of local atoms (atomCapacity), or the specified number
      * of ghosts (ghostCapacity), and allocate send and receive 
      * buffers large enough to hold the larger value.
      *
      * \param atomCapacity max expected number of local atoms.
      * \param ghostCapacity max expected number of ghost atoms.
      */
      void allocate(int atomCapacity, int ghostCapacity);

      /// \name Send Buffer Management
      ///@{
      
      /**
      * Clear the send buffer.
      */
      void clearSendBuffer();

      /**
      * Initialize a data block.
      *
      * Sets sendSize() to zero and sets the sendType.
      *
      * \param sendType BlockDataType value for type of data to be sent
      */
      void beginSendBlock(int sendType);

      /**
      * Function template for packing one variable into the send buffer.
      *
      * This template is used to implement the pack functions provided
      * by the DdMd::Atom and DdMd::Group classes. It is designed to 
      * pack a single primitive C variable into the send buffer. It will 
      * also work on any plain old data type T for which assignment (=) 
      * does a straight bitwise copy.
      *
      * \param data variable to be packed into the send buffer.
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
      * This method writes the "prefix" section at the beginning of a data
      * block. The descriptor specifies the length of the block, the type 
      * of data, and whether the block is "complete". 
      *
      * A block should be marked as incomplete iff the required data of 
      * the relevant sendtype did not fit into the buffer. This tells the
      * receiving processor to expect one or more other buffers containing 
      * the remaining data of the same type type.
      *
      * \param isComplete false if data block is incomplete, true otherwise
      */
      void endSendBlock(bool isComplete = true);

      ///@}
      /// \name Receive Buffer Management
      ///@{
      
      /**
      * Begin to receive a block from the recv buffer.
      *
      * This method reads the descriptor section at the beginning of a
      * block, and sets the receive cursor to be ready to unpack the 
      * first item.
      *
      * \return false if this block is complete, false otherwise.
      */
      bool beginRecvBlock();

      /**
      * Function template unpacking one variable from the receive buffer.
      *
      * This template is used to implement the unpack functions provided
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
      * Decrement recvSize counter after completely unpacking an item.
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

      ///@}
      #ifdef UTIL_MPI
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
      #endif
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
      * Call this on master, after calling computeStatistics on all processors.
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
      * Number of unread items left in current recv block.
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
