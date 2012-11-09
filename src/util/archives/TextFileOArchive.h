#ifndef UTIL_TEXT_FILE_O_ARCHIVE_H
#define UTIL_TEXT_FILE_O_ARCHIVE_H

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
   * Saving archive for character based ostream.
   *
   * \ingroup Archive_Module
   */
   class TextFileOArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      TextFileOArchive();

      /**
      * Destructor.
      */
      virtual ~TextFileOArchive();

      /**
      * Set the stream.
      *
      * \param out output stream to which to write.
      */
      void setStream(std::ostream& out);

      /**
      * Serialize one object to this archive.
      */
      template <typename T>
      TextFileOArchive& operator & (T& data);

      /**
      * Serialize one object.
      */
      template <typename T>
      TextFileOArchive& operator << (T& data);

      /**
      * Serialize a single T object to file as formatted ascii.
      *
      * \param data object to be written to file
      */
      template <typename T> 
      void pack(const T& data);

      /**
      * Serialize an array of objects to file.
      *
      * \param array C-array of T objects
      * \param n     number of elements
      */
      template <typename T> 
      void pack(const T* array, int n);

   private:

      /// Pointer to output stream file.
      std::ostream* ostreamPtr_;

      /// Archive version id.
      unsigned int  version_;

   };

   // Inline methods

   inline bool TextFileOArchive::is_saving()
   {  return true; }

   inline bool TextFileOArchive::is_loading()
   {  return false; }

   // Inline non-static method templates.

   /*
   * Serialize one object.
   */
   template <typename T>
   inline TextFileOArchive& TextFileOArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Serialize one object to this archive.
   */
   template <typename T>
   inline TextFileOArchive& TextFileOArchive::operator << (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Method templates

   /*
   * Serialize a single object of type T.
   */
   template <typename T>
   inline void TextFileOArchive::pack(const T& data)
   {  *ostreamPtr_ << data << std::endl; }

   /*
   * Serialize a single object of type double.
   */
   template <>
   inline void TextFileOArchive::pack(const double& data)
   {  
      ostreamPtr_->setf(std::ios::scientific);
      ostreamPtr_->width(25);
      ostreamPtr_->precision(17);
      *ostreamPtr_ << data << std::endl; 
   }

   /*
   * Serialize a C-array of objects of type T.
   */
   template <typename T>
   inline void TextFileOArchive::pack(const T* array, int n)
   {
      for (int i=0; i < n; ++i) {
        *ostreamPtr_ << array[i] << "  ";
      }
      *ostreamPtr_ << std::endl;
   }

   /*
   * Serialize a C-array of double values.
   */
   template <>
   inline void TextFileOArchive::pack(const double* array, int n)
   {
      ostreamPtr_->setf(std::ios::scientific);
      ostreamPtr_->precision(17);
      for (int i=0; i < n; ++i) {
        ostreamPtr_->width(25);
        *ostreamPtr_ << array[i] << "  ";
      }
      *ostreamPtr_ << std::endl;
   }

   // Explicit serialize functions for primitive types

   /*
   * Serialize a bool to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize a char to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize an unsigned int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize an int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize an unsigned long int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize a long int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize a float to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize an double to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   // Explicit serialize functions for primitive types

   /*
   * Serialize a std::complex<float> to a TextFileOArchive.
   */
   template <>
   inline 
   void serialize(TextFileOArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize a std::complex<double> to a TextFileOArchive.
   */
   template <>
   inline 
   void serialize(TextFileOArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /*
   * Serialize a std::string to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, std::string& data, 
                         const unsigned int version)
   {
      int size = data.size();
      ar.pack(size);
      ar.pack(data);
   }

   // Explicit serialize functions for namespace Util

   /*
   * Serialize a Util::Vector to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.pack(data); } 

   /*
   * Serialize a Util::IntVector to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.pack(data); }

}
#endif
