#ifndef UTIL_BINARY_FILE_I_ARCHIVE_H
#define UTIL_BINARY_FILE_I_ARCHIVE_H

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
   * Saving archive for binary istream.
   *
   * \ingroup Archive_Module
   */
   class BinaryFileIArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      BinaryFileIArchive();

      /**
      * Destructor.
      */
      virtual ~BinaryFileIArchive();

      /**
      * Set the stream.
      *
      * \param out output stream to which to write.
      */
      void setStream(std::istream& out);

      /**
      * Read one object.
      */
      template <typename T>
      BinaryFileIArchive& operator & (T& data);

      /**
      * Read one object.
      */
      template <typename T>
      BinaryFileIArchive& operator >> (T& data);

      template <typename T> 
      void unpack(T& data);

      template <typename T> 
      void unpack(T* array, int n);

   private:

      /// Pointer to output stream file.
      std::istream* istreamPtr_;

      /// Archive version id.
      unsigned int  version_;

   };

   // Inline methods

   inline bool BinaryFileIArchive::is_saving()
   {  return false; }

   inline bool BinaryFileIArchive::is_loading()
   {  return true; }

   // Inline non-static method templates

   /*
   * Read one object.
   */
   template <typename T>
   inline BinaryFileIArchive& BinaryFileIArchive::operator & (T& data)
   {  
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Read one object.
   */
   template <typename T>
   inline BinaryFileIArchive& BinaryFileIArchive::operator >> (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Method templates

   /*
   * Load a single object of type T.
   */
   template <typename T>
   inline void BinaryFileIArchive::unpack(T& data)
   {  istreamPtr_->read( (char*)(&data), sizeof(T) ); }

   /*
   * Load a C-array of objects of type T.
   */
   template <typename T>
   inline void BinaryFileIArchive::unpack(T* array, int n)
   {
      for (int i=0; i < n; ++i) {
         istreamPtr_->read( (char*)(&array[i]), sizeof(T));
      }
   }

   // Explicit serialize functions for primitive types

   /*
   * Load a bool from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a char from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned int from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an int from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned long int from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a long int from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a float from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a double from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   // Explicit serialize functions for std library types

   /*
   * Load a std::complex<float> from a BinaryFileIArchive.
   */
   template <>
   inline 
   void serialize(BinaryFileIArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::complex<double> from a BinaryFileIArchive.
   */
   template <>
   inline 
   void serialize(BinaryFileIArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::string from a BinaryFileIArchive.
   */
   template <>
   void serialize(BinaryFileIArchive& ar, std::string& data, 
                  const unsigned int version);

   // Explicit serialize functions for namespace Util types

   /*
   * Load a Util::Vector from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.unpack(data); } 

   /*
   * Load a Util::IntVector from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

}
#endif
