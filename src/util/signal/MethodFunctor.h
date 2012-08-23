#ifndef METHOD_FUNCTOR_H
#define METHOD_FUNCTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "IFunctor.h"

namespace Util
{

   /**
   * Functor that wraps a one-argument class member function.
   *
   * The constructor to MethodFunctor<T> takes pointers to an invoking instance
   * of class Object and a member function (method) that takes one T argument.
   * The operator (const T&) invokes that method on that object.
   */
   template <class Object, typename T=void>
   class MethodFunctor : public IFunctor<T>
   {
   public:

      typedef void (Object::*Method1Ptr)(const T&);

      /**
      * Constructor.
      *
      * \param objectPtr pointer to invoking object
      * \param methodPtr pointer to member function
      */
      MethodFunctor(Object* objectPtr, Method1Ptr methodPtr) 
       : objectPtr_(objectPtr),
         methodPtr_(methodPtr)
      {}

      virtual void operator () (const T& t)
      {  (objectPtr_->*methodPtr_)(t); }

   private:

      Object*     objectPtr_;
      Method1Ptr  methodPtr_;

   };

   /**
   * Functor that wraps a class member function with no arguments.
   */
   template <class Object>
   class MethodFunctor<Object, void> : public IFunctor<void>
   {
   public:

      typedef void (Object::*Method0Ptr)();

      /**
      * Constructor.
      *
      * \param objectPtr pointer to invoking object
      * \param methodPtr pointer to member function
      */
      MethodFunctor(Object* objectPtr, Method0Ptr methodPtr) 
       : objectPtr_(objectPtr),
         methodPtr_(methodPtr)
      {}

      virtual void operator () ()
      {  (objectPtr_->*methodPtr_)(); }

   private:

      Object*    objectPtr_;
      Method0Ptr methodPtr_;

   };

}
#endif 
