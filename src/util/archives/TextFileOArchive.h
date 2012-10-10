#ifndef TEXT_FILE_O_ARCHIVE_H
#define TEXT_FILE_O_ARCHIVE_H

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
#include <fstream>

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
      * Constructor.
      *
      * \param filename name of file to open for reading.
      */
      TextFileOArchive(std::string filename);

      /**
      * Destructor.
      */
      virtual ~TextFileOArchive();

      /**
      * Get the underlying ifstream by reference.
      */
      std::ofstream& file();

      /**
      * Write one object.
      */
      template <typename T>
      TextFileOArchive& operator & (T& data);

      /**
      * Write one object.
      */
      template <typename T>
      TextFileOArchive& operator << (T& data);

      /**
      * Pack one object.
      */
      template <typename T> 
      void pack(const T& data);

      /**
      * Pack a C array.
      *
      * \param array address of a C array (or first element)
      * \param n     number of elements
      */
      template <typename T> 
      void pack(const T* array, int n);

      /**
      * Pack a 2D C array.
      *
      * \param array address of first row (or element) of 2D array
      * \param m     number of rows
      * \param n     number of columns
      */
      template <typename T> 
      void pack(const T* array, int m, int n);

   private:

      /// Pointer to output stream file.
      std::ofstream* filePtr_;

      /// Archive version id.
      unsigned int  version_;

   };

   // Inline methods

   inline bool TextFileOArchive::is_saving()
   {  return true; }

   inline bool TextFileOArchive::is_loading()
   {  return false; }

   // Inline non-static methods

   /*
   * Write one object.
   */
   template <typename T>
   inline TextFileOArchive& TextFileOArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Write one object.
   */
   template <typename T>
   inline TextFileOArchive& TextFileOArchive::operator << (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Method templates

   /*
   * Pack a single object of type T.
   */
   template <typename T>
   inline void TextFileOArchive::pack(const T& data)
   {  *filePtr_ << data << std::endl; }

   /*
   * Pack a single object of type double.
   */
   template <>
   inline void TextFileOArchive::pack(const double& data)
   {  
      filePtr_->setf(std::ios::scientific);
      filePtr_->width(25);
      filePtr_->precision(17);
      *filePtr_ << data << std::endl; 
   }

   /*
   * Write a C-array of objects of type T.
   */
   template <typename T>
   inline void TextFileOArchive::pack(const T* array, int n)
   {
      for (int i=0; i < n; ++i) {
        *filePtr_ << array[i] << "  ";
      }
      *filePtr_ << std::endl;
   }

   /*
   * Pack a 2D C-array of objects of type T.
   */
   template <typename T>
   inline void TextFileOArchive::pack(const T* array, int m, int n)
   {
      int i, j;
      for (i=0; i < m; ++i) {
         for (j=0; j < n; ++j) {
            *filePtr_ << array[i*n + j] << "  ";
         }
         *filePtr_ << std::endl;
      }
   }

   // Explicit serialize functions for primitive types

   /**
   * Save a bool to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save a char to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save an unsigned int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save an int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save an unsigned long int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save a long int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save a float to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save an double to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   // Explicit serialize functions for primitive types

   /**
   * Save a std::complex<float> to a TextFileOArchive.
   */
   template <>
   inline 
   void serialize(TextFileOArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save a std::complex<double> to a TextFileOArchive.
   */
   template <>
   inline 
   void serialize(TextFileOArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /**
   * Save a std::string to a TextFileOArchive.
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

   /**
   * Save a Util::Vector to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.pack(data); } 

   /**
   * Save a Util::IntVector to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.pack(data); }

}
#endif
