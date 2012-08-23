#ifndef SIGNAL_H
#define SIGNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <list>

namespace Util
{

   using std::list;

   // Functions with one argument

   /**
   * Functor that wraps member function with one argument (abstract).
   */
   class Signal1 {

   public: 

      /**
      * Call a specific member function.
      */
      virtual void operator () (T1&) = 0; 

   }

   /**
   * Implementation of functor for member function with 1 argument.
   */
   template <typename Object, typename T1>
   class Signal1Impl
   {
   public:

      /**
      * Constructor.
      *
      * \param objectPtr pointer to invoking object
      * \param memberPtr pointer to member function
      */
      Single1Impl(Object* objectPtr, Object::*memberPtr(T1& )) 
       : objectPtr_(objectPtr),
         memberPtr_(memberPtr)
      {}

      virtual void notify(T1& t1)
      {  objectPtr_->*memberPtr_(t1); }

   private:

      Object*    objectPtr_;
      Object::*  memberPtr_;

   };

}
#endif 
