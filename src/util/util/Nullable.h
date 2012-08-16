#ifndef NULLABLE_H
#define NULLABLE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * Template for a value that can be declared null.
   *
   * \ingroup Param_Module
   */
   template <class T>
   class Nullable
   {

   public:

      /**
      * Default constructor.
      */
      Nullable()
       : value_(),
         isNull_(true)
      {}

      /**
      * Copy constructor.
      *
      * \param other Nullable object being copied.
      */
      Nullable(const Nullable<T>& other)
       : value_(other.value_),
         isNull_(other.isNull_)
      {}

      /**
      * Construct from T value (explciit).
      *
      * \param value value of wrapped object.
      */
      explicit Nullable(const T& value)
       : value_(value),
         isNull_(false)
      {}

      /**
      * Assignment.
      */ 
      Nullable<T>& operator = (const Nullable<T>& other)
      {
         if (this != &other) {
            value_  = other.value_;
            isNull_ = other.isNull_;
         }
         return *this;
      }

      /**
      * Assignment from value.
      */ 
      Nullable<T>& operator = (const T& value)
      {
         value_  = value;
         isNull_ = false;
         return *this;
      }

      /**
      * Mark value as null (unknown).
      */
      void setNull()
      {  isNull_ = true; }

      /**
      * Is this object null (is the value unknown)?
      */
      bool isNull()
      {  return isNull_; }

      /**
      * Return value (if not null)
      */
      T& value()
      {
         if (isNull_) {
            UTIL_THROW("Attempt to return null value.");
         }  
         return value_; 
      }

      #if 0
      /**
      * Conversion to T value, if not null.
      *
      * Throws an exception if this object is null.
      */
      operator T ()
      {
         if (isNull_) {
            UTIL_THROW("Attempted conversion of null value");
         }
         return value_;
      }
      #endif

   private:

      /// Value of associated object.
      T   value_;

      /// Set true if value is unknown (null).
      bool isNull_;

   };

} 
#endif
