#ifndef PACKED_DATA_H
#define PACKED_DATA_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * Memory block for packed heterogeneous binary data.
   *
   * A PackedData object encapsulates a unsigned char array 
   * that may be used to store packed, heterogenous binary 
   * data. A PackedData object must be allocated before it 
   * can be used.
   *
   * To pack data, first call beginPacking(), followed by a
   * series of calls of the pack methods. During packing, the
   * cursor() method returns the address of the next available 
   * unwritten byte.  The overloaded pack method templates write 
   * either a single object or an array of objects to memory at 
   * the current cursor location, and then updates the cursor. 
   * Packing is terminated by calling beginUnpacking(), which
   * sets the cursor to the beginning in preparation for
   * unpacking. 
   *
   * To unpack data, call beginUnpacking(), followed by a series
   * of calls of the unpack methods.  During unpacking, cursor()
   * points to the first unread byte, and endPacked() points to
   * end of the packed block (one byte after the last). Unpacking 
   * is terminated by calling beginPacking(), which resets the 
   * cursor to the beginning in preparation for re-packing. 
   *
   * The pack and unpack methods may not be intermingled. One or 
   * more calls to pack must be preceded by a call of beginPacking(), 
   * and one or more calls to unpack must be preceded by a call to 
   * beginUnpacking().
   *
   * Example:
   * \code
   *
   * PackedData buffer;
   * int capacity = 128;
   * buffer.allocate(capacity);
   *
   * const int n  = 2;
   *
   * // Variables to be packed int buffer
   * int i = 3;
   * double d = 5.0;
   * double a[n]  a;
   * a[0] = 9.0;
   * a[1] = 8.0;
   *
   * buffer.beginPacking();
   * buffer.pack(i);
   * buffer.pack(d);
   * buffer.pack(a, n); 
   *
   * // Variables to be unpacked from buffer
   * int j;
   * double e;
   * double[n] b;
   *
   * buffer.beginUnpacking();
   * buffer.unpack(j);
   * buffer.unpack(e);
   * buffer.unpack(b, n); 
   *
   * \endcode
   */
   class PackedData
   {

   public:

      /**
      * Define a "Byte" type.
      */
      typedef char Byte;

      /**
      * Constructor.
      */
      PackedData();

      /**
      * Destructor.
      */
      ~PackedData();

      /**
      * Allocate memory block.
      *
      * \param capacity sizeof of block, in Bytes.
      */
      void allocate(int capacity);

      /**
      * Begin packing, set cursor to the beginning.
      */
      void beginPacking();

      /**
      * Begin unpacking, set cursor to the beginning.
      */
      void beginUnpacking();

      /**
      * Pack one object.
      *
      * Upon return, the cursor points to next unwritten Byte.
      *
      * \param data variable to be packed.
      */
      template <typename T>
      void pack(const T& data);

      /**
      * Pack a C array of objects of type T.
      *
      * Upon return, the cursor points to next unwritten Byte.
      *
      * \param array pointer to first element in array.
      * \param n     number of elements in array.
      */
      template <typename T>
      void pack(const T* array, int n);

      /**
      * Unpack one object.
      *
      * Upon return, the cursor points to next unread Byte.
      *
      * \param data object to be unpacked
      */
      template <typename T>
      void unpack(T& data);

      /**
      * Unpack a C array and update cursor .
      *
      * Upon return, the cursor points to next unread Byte.
      *
      * \param array pointer to first element in array.
      * \param n     number of elements in array.
      */
      template <typename T>
      void unpack(T* array, int n);

      #if 0
      #ifdef UTIL_MPI
      /**
      * Send packed data via MPI.
      *
      * \param comm  MPI communicator
      * \param dest  rank of processor to which data is sent
      */
      void send(MPI::Intracomm& comm, int dest);

      /**
      * Receive packed data via MPI.
      *
      * \param comm   MPI communicator
      * \param source rank of processor from which data is sent.
      */
      void recv(MPI::Intracomm& comm, int source);
      #endif
      #endif

      /**
      * Return pointer to beginning of block.
      */
      Byte* begin() const;

      /**
      * Return pointer to end of allocated block (one Byte past the last).
      */
      Byte* endAllocated() const;

      /**
      * Return pointer to end of packed block (one Byte past the last).
      */
      Byte* endPacked() const;

      /**
      * Return pointer to cursor position (cursor).
      */
      Byte* cursor() const;

      /**
      * Has memory been allocated?
      */
      bool isAllocated() const;

      /**
      * Return capacity in Bytes.
      */
      int capacity() const;

   private:

      /**
      * Enumeration for current status.
      */
      enum Status {NONE, PACKING, UNPACKING};

      /// Pointer to first element in block. 
      Byte* begin_;

      /// Pointer to one Byte past last in allocated block.
      Byte* endAllocated_;

      /// Pointer to one Byte past last byte in packed block.
      Byte* endPacked_;

      /// Current element (read/write cursor).
      Byte* cursor_;

      /// Allocated size of send and recv buffers, in Bytes.
      int capacity_;

      /// Is the buffer unused, packing, or unpacking?
      Status status_;

   };

   // Inline methods

   /*
   * Has this data block been allocated?
   */
   inline bool PackedData::isAllocated() const
   {  return (bool) begin_; }

   /*
   * Return capacity in Bytes.
   */
   inline int PackedData::capacity() const
   {  return capacity_; }

   /*
   * Return pointer to beginning of block.
   */
   inline PackedData::Byte* PackedData::begin() const
   {  return begin_; }

   /*
   * Return end of allocated block (one Byte past the last).
   */
   inline PackedData::Byte* PackedData::endAllocated() const
   {  return endAllocated_; }

   /*
   * Return end of packed block (one Byte past the last).
   */
   inline PackedData::Byte* PackedData::endPacked() const
   {  return endPacked_; }

   /*
   * Return pointer to cursor position.
   */
   inline PackedData::Byte* PackedData::cursor() const
   {  return cursor_; }

   // Template methods

   /*
   * Pack an object of type T into send block.
   */
   template <typename T>
   void PackedData::pack(const T& data)
   {
      if (status_ != PACKING) {
         UTIL_THROW("Status is not PACKING");
      }
      if (cursor_ + sizeof(T) > endAllocated_) {
         UTIL_THROW("Attempted write past end of send block");
      }
      T* ptr = (T *)cursor_;
      *ptr = data;
      ++ptr;
      cursor_ = (Byte *)ptr;
   }

   /*
   * Pack a C-array of objects of type T.
   */
   template <typename T>
   void PackedData::pack(const T* array, int n)
   {
      if (status_ != PACKING) {
         UTIL_THROW("Status is not PACKING");
      }
      if (cursor_ + n*sizeof(T) > endAllocated_) {
         UTIL_THROW("Attempted write past end of send block");
      }
      T* ptr = (T *)cursor_;
      for (int i=0; i < n; ++i) {
         *ptr = array[i];
         ++ptr;
      }
      cursor_ = (Byte *)ptr;
   }

   /*
   * Unpack one object of type T.
   */
   template <typename T>
   void PackedData::unpack(T& data)
   {
      if (status_ != UNPACKING) {
         UTIL_THROW("Status is not UNPACKING");
      }
      if (cursor_ + sizeof(data) > endPacked_) {
         UTIL_THROW("Attempted read past end of packed block");
      }
      T* ptr = (T *)cursor_;
      data = *ptr;
      ++ptr;
      cursor_ = (Byte *)ptr;
   }

   /*
   * Unpack a C-array of objects of type T.
   */
   template <typename T>
   void PackedData::unpack(T* array, int n)
   {
      if (status_ != UNPACKING) {
         UTIL_THROW("Status is not UNPACKING");
      }
      if (cursor_ + n*sizeof(T) > endPacked_) {
         UTIL_THROW("Attempted write past end of send block");
      }
      T* ptr = (T *)cursor_;
      for (int i=0; i < n; ++i) {
         array[i] = *ptr;
         ++ptr;
      }
      cursor_ = (Byte *)ptr;
   }

}
#endif
