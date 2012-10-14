#ifndef UTIL_BINARY_FILE_O_ARCHIVE_H
#define UTIL_BINARY_FILE_O_ARCHIVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Byte.h"
#include "serialize.h"

#include <util/space/Vector.h>
#include <util/space/IntVector.h>

#include <complex>
#include <string>
#include <iostream>

namespace Util
{

   /**
   * Saving / output archive for binary ostream.
   *
   * \ingroup Archive_Module
   */
   class BinaryFileOArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      BinaryFileOArchive();

      /**
      * Destructor.
      */
      virtual ~BinaryFileOArchive();

      /**
      * Set the stream.
      *
      * \param out output stream to which to write.
      */
      void setStream(std::ostream& out);

      /**
      * Write one object.
      */
      template <typename T>
      BinaryFileOArchive& operator & (T& data);

      /**
      * Write one object.
      */
      template <typename T>
      BinaryFileOArchive& operator << (T& data);

      template <typename T> 
      void pack(const T& data);

      template <typename T> 
      void pack(const T* array, int n);

   private:

      /// Pointer to output stream file.
      std::ostream* ostreamPtr_;

      /// Archive version id.
      unsigned int  version_;

   };

   // Inline methods

   inline bool BinaryFileOArchive::is_saving()
   {  return true; }

   inline bool BinaryFileOArchive::is_loading()
   {  return false; }

   // Inline non-static methods

   /*
   * Write one object.
   */
   template <typename T>
   inline BinaryFileOArchive& BinaryFileOArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Write one object.
   */
   template <typename T>
   inline BinaryFileOArchive& BinaryFileOArchive::operator << (T& data)
   {  
      serialize(*this, data, version_);
      return *this;
   }

   // Method templates

   /*
   * Bitwise pack a single object of type T.
   */
   template <typename T>
   inline void BinaryFileOArchive::pack(const T& data)
   {  ostreamPtr_->write( (char*)(&data), sizeof(T) ); }

   /*
   * Bitwise pack a C-array of objects of type T.
   */
   template <typename T>
   inline void BinaryFileOArchive::pack(const T* array, int n)
   {
      for (int i=0; i < n; ++i) {
         ostreamPtr_->write( (char*)(&array[i]), sizeof(T));
      }
   }

   // Explicit serialize functions for primitive types

   /*
   * Save a bool to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a char to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an unsigned int to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an int to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an unsigned long int to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a long int to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a float to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an double to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   // Explicit serialize functions for std library types

   /*
   * Save a std::complex<float> to a BinaryFileOArchive.
   */
   template <>
   inline 
   void serialize(BinaryFileOArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a std::complex<double> to a BinaryFileOArchive.
   */
   template <>
   inline 
   void serialize(BinaryFileOArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a std::string to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, std::string& data, 
                         const unsigned int version)
   {
      size_t size = data.size() + 1; // the +1 is for the NULL
      ar.pack(size);
      const char* temp = data.c_str();
      ar.pack(temp, size);
   }

   // Explicit serialize functions for namespace Util

   /*
   * Save a Util::Vector to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.pack(data); } 

   /*
   * Save a Util::IntVector to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.pack(data); }

}
#endif
