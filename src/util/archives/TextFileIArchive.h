#ifndef UTIL_TEXT_FILE_I_ARCHIVE_H
#define UTIL_TEXT_FILE_I_ARCHIVE_H

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
   * Loading archive for text istream.
   *
   * \ingroup Archive_Module
   */
   class TextFileIArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      TextFileIArchive();

      /**
      * Destructor.
      */
      virtual ~TextFileIArchive();

      /**
      * Set the stream.
      *
      * \param out output stream to which to write.
      */
      void setStream(std::istream& out);

      /**
      * Load one object.
      */
      template <typename T>
      TextFileIArchive& operator & (T& data);

      /**
      * Load one object.
      */
      template <typename T>
      TextFileIArchive& operator >> (T& data);

      /**
      * Serialize a single T object.
      *
      * \param data object to be serialized.
      */
      template <typename T> 
      void unpack(T& data);

      /**
      * Serialize a C-array of T objects.
      *
      * \param array pointer to array of T objecs.
      * \param n     number of elements in array
      */
      template <typename T> 
      void unpack(T* array, int n);

   private:

      /// Pointer to output stream file.
      std::istream* istreamPtr_;

      /// Archive version id.
      unsigned int  version_;

   };

   // Inline methods

   inline bool TextFileIArchive::is_saving()
   {  return false; }

   inline bool TextFileIArchive::is_loading()
   {  return true; }

   // Inline non-static methods

   /*
   * Load one object.
   */
   template <typename T>
   inline TextFileIArchive& TextFileIArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Load one object.
   */
   template <typename T>
   inline TextFileIArchive& TextFileIArchive::operator >> (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Method templates

   /*
   * Load a single object of type T.
   */
   template <typename T>
   inline void TextFileIArchive::unpack(T& data)
   {  *istreamPtr_ >> data; }

   /*
   * Load a C array.
   */
   template <typename T>
   inline void TextFileIArchive::unpack(T* array, int n)
   {
      for (int i=0; i < n; ++i) {
         *istreamPtr_ >> array[i];
      }
   }

   /*
   * Unpack a single char.
   */
   template <>
   inline void TextFileIArchive::unpack(char& data)
   {   istreamPtr_->get(data); }

   /*
   * Unpack a C-array of characters.
   */
   template <>
   inline void TextFileIArchive::unpack(char* array, int n)
   {
      istreamPtr_->get(array, n+1,'\0');
   }

   // Explicit serialize functions for primitive types

   /*
   * Load a bool from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a char from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned long int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a long int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a float from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a double from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   // Explicit serialize functions for primitive types

   /*
   * Load a std::complex<float> from a TextFileIArchive.
   */
   template <>
   inline 
   void serialize(TextFileIArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::complex<double> from a TextFileIArchive.
   */
   template <>
   inline 
   void serialize(TextFileIArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::string from a TextFileIArchive.
   */
   template <>
   void serialize(TextFileIArchive& ar, std::string& data, 
                  const unsigned int version);

   // Explicit serialize functions for namespace Util types

   /*
   * Load a Util::Vector from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.unpack(data); } 

   /*
   * Load a Util::IntVector from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

}
#endif
