#ifndef BUFFER_H
#define BUFFER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/
#include <util/param/ParamComposite.h>  // base class
#include <ddMd/chemistry/Bond.h>        // typedef
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   class Atom;
   template <int N> class Group;

   /**
   * Buffer for sending local or ghost atoms.
   */
   class Buffer: public ParamComposite 
   {

   public:

      enum BlockDataType {NONE, ATOM, GHOST, GROUP};

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
      void readParam(std::istream& in);

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

      #ifdef UTIL_MPI
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
      * Finalize a block in the send buffer.
      *
      * \param isComplete false if data block is incomplete, true otherwise.
      */
      void endSendBlock(bool isComplete = true);

      /**
      * Begin to receive a block from the recv buffer.
      *
      * \return false if an incomplete block was received, true otherwise.
      */
      bool beginRecvBlock();

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

      /**
      * Pack an Atom for exchange of ownership.
      *
      * Appends an atom to the send buffer, and increments sendSize().
      *
      * Throws an Exception if sendType != ATOM.
      *
      * \param atom Atom object to be sent.
      */
      void packAtom(Atom& atom);

      /**
      * Receive ownership of an Atom.
      *
      * Unpacks an atom from recv buffer into atom, and decrements 
      * recvSize().
      *
      * \param atom Atom object into which data is received.
      */
      void unpackAtom(Atom& atom);

      /**
      * Pack ghost Atom into send buffer.
      *
      * Copies required data into send buffer, and increments sendSize().
      *
      * \param atom ghost Atom to be sent.
      */
      void packGhost(Atom& atom);

      /**
      * Unpack a ghost Atom from recv buffer.
      *
      * Unpacks required data from recv buffer into atom, decrements
      * recvSize().
      *
      * \param atom ghost Atom object.
      */
      void unpackGhost(Atom& atom);

      /**
      * Pack a Group for sending.
      *
      * \param group Group<N> object to be sent.
      */
      template <int N>
      void packGroup(const Group<N>& group);

      /**
      * Receive the next group in the recv buffer.
      *
      * \param group Group<N> object into which data is received.
      */
      template <int N>
      void unpackGroup(Group<N>& group);

      /**
      * Discard the next group in the recv buffer.
      */
      template <int N>
      void discardGroup();

      #if 0
      /**
      * Pack a Bond into send buffer.
      *
      * Copies required data from atom into send buffer, and
      * increments sendSize() and send buffer pointer.
      *
      * Throws an Exception if sendType != BOND.
      *
      * \param bond Bond to be sent.
      */
      void packBond(Bond& bond);

      /**
      * Unpack a Bond from recv buffer.
      *
      * Copies required data from recv buffer into bond, and
      * decrements recvSize() and recv buffer pointer.
      *
      * \param atom ghost Atom object.
      */
      void unpackBond(Bond& bond);
      #endif

      /**
      * Number of items written to current send block.
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

      #endif

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

   private:

      #ifdef UTIL_MPI
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

      /// Allocated size of send and recv buffers, in bytes.
      int bufferCapacity_;

      /// Number of atoms or ghosts currently in send buffer.
      int sendSize_;

      /// Number of unread atoms or ghosts currently in receive buffer.
      int recvSize_;

      /// Number of atoms in group (or 0 if not a Group).
      int sendGroupSize_;

      /// Type of atom being sent = NONE, ATOM, GHOST, BOND, ANGLE. DIHEDRAL
      BlockDataType sendType_;

      /// Type of atom being received (BlockDataType cast to int)
      int recvType_;

      /// Number of atoms in group (or 0 if not a Group).
      int recvGroupSize_;
      #endif

      /// Maximum number of local atoms in buffer.
      int atomCapacity_;

      /// Maximum number of ghost atoms in buffer.
      int ghostCapacity_;

      /// Has this buffer been initialized ?
      bool isInitialized_;

      /// Number of bytes required per local atom.
      static const int localAtomSize_ = 60;

      /// Number of bytes required per ghost atom.
      static const int ghostAtomSize_ = 36;

      #ifdef UTIL_MPI
      /**
      * Template method for packing a variable
      *
      * \param data variable to be packed.
      */
      template <typename T>
      void pack(const T& data);

      /**
      * Template method for unpacking a variable
      *
      * \param data variable to be unpacked.
      */
      template <typename T>
      void unpack(T& data);

      /*
      * Allocate send and recv buffers, using preset capacities.
      */
      void allocate();
      #endif

   };

   #if UTIL_MPI
   /*
   * Pack an object of type T into send buffer.
   */
   template <typename T>
   void Buffer::pack(const T& data)
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
   void Buffer::unpack(T& data)
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
   * Pack a Group.
   */
   template <int N>
   void Buffer::packGroup(const Group<N>& group)
   {
      // Preconditions
      if (sendType_ != GROUP) {
         UTIL_THROW("SendType is not GROUP");
      }
      if (sendGroupSize_ != N) {
         UTIL_THROW("Wrong sendGroupSize");
      }
      pack<int>(group.typeId());
      pack<int>(group.id());
      for (int j = 0; j < N; ++j) {
         pack<int>(group.atomId(j));
      }

      //Increment number of groups in send buffer by 1
      ++sendSize_;

   }

   /*
   * Unpack the next Group in the recieve buffer.
   */
   template <int N>
   void Buffer::unpackGroup(Group<N>& group)
   {
      // Preconditions
      if (recvType_ != GROUP) {
         UTIL_THROW("SendType is not GROUP");
      }
      if (recvGroupSize_ != N) {
         UTIL_THROW("Wrong recvGroupSize");
      }
      int i, j;
      unpack(i);
      group.setTypeId(i);
      unpack(i);
      group.setId(i);
      for (j = 0; j < N; ++j) {
         unpack(i);
         group.setAtomId(j, i);
         group.clearAtomPtr(j);
      }
      group.setPostMark(false);

      // Decrement number of groups in recv buffer by 1
      recvSize_--;

   }

   /*
   * Discard the next Group in the receive buffer.
   */
   template <int N>
   void Buffer::discardGroup()
   {
      // Preconditions
      if (recvType_ != GROUP) {
         UTIL_THROW("SendType is not GROUP");
      }
      if (recvGroupSize_ != N) {
         UTIL_THROW("Wrong recvGroupSize");
      }
      int i, j;
      unpack(i);
      unpack(i);
      for (j = 0; j < N; ++j) {
         unpack(i);
      }

      // Decrement number of groups in recv buffer by 1
      recvSize_--;

   }
   #endif

}
#endif
