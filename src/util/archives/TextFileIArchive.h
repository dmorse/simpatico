#ifndef TEXT_FILE_I_ARCHIVE_H
#define TEXT_FILE_I_ARCHIVE_H

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
      * Write one object.
      */
      template <typename T>
      TextFileIArchive& operator & (T& data);

      /**
      * Write one object.
      */
      template <typename T>
      TextFileIArchive& operator >> (T& data);

      /**
      * Unpack an item of type T.
      */
      template <typename T> 
      void unpack(T& data);

      /**
      * Unpack a C array.
      *
      * \param array C array
      * \param n     number of elements
      */
      template <typename T> 
      void unpack(T* array, int n);

      /**
      * Unpack a 2D C array.
      *
      * \param array pointer to first element
      * \param m     number of rows
      * \param n     number of columns
      */
      template <typename T> 
      void unpack(T* array, int m, int n);

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
   * Write one object.
   */
   template <typename T>
   inline TextFileIArchive& TextFileIArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Write one object.
   */
   template <typename T>
   inline TextFileIArchive& TextFileIArchive::operator >> (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Method templates

   /*
   * Bitwise unpack a single object of type T.
   */
   template <typename T>
   inline void TextFileIArchive::unpack(T& data)
   {  *istreamPtr_ >> data; }

   /*
   * Unpack a C array.
   */
   template <typename T>
   inline void TextFileIArchive::unpack(T* array, int n)
   {
      for (int i=0; i < n; ++i) {
         *istreamPtr_ >> array[i];
      }
   }

   /*
   * Bitwise pack a 2D C-array of objects of type T.
   */
   template <typename T>
   void TextFileIArchive::unpack(T* array, int m, int n)
   {
      int i, j;
      for (i = 0; i < m; ++i) {
         for (j = 0; j < n; ++j) {
            *istreamPtr_ >> array[i*n + j];
         }
      }
   }

   /*
   * Read a single char.
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

   /**
   * Load a bool from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Load a char from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Save an unsigned int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Save an int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Load an unsigned long int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Load a long int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Load a float from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Save an double from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   // Explicit serialize functions for primitive types

   /**
   * Load a std::complex<float> from a TextFileIArchive.
   */
   template <>
   inline 
   void serialize(TextFileIArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Load a std::complex<double> from a TextFileIArchive.
   */
   template <>
   inline 
   void serialize(TextFileIArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /**
   * Load a std::string from a TextFileIArchive.
   */
   template <>
   void serialize(TextFileIArchive& ar, std::string& data, 
                  const unsigned int version);

   // Explicit serialize functions for namespace Util

   /**
   * Load a Util::Vector from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.unpack(data); } 

   /**
   * Load a Util::IntVector from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

}
#endif
