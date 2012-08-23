#ifndef SIGNAL_H
#define SIGNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "IFunctor.h"
#include "MethodFunctor.h"

#include <list>

namespace Util
{

   using std::list;

   // Signal with one argument.

   /**
   * Notifier (or subject) in the Observer design pattern.
   *
   * A Signal manages a list of registered functor objects, and provides
   * a notify() method that calls them all with the same argument.
   * 
   * \ingroup Util_Module
   */
   template <typename T=void>
   class Signal
   {
   
   public:

      // Compiler default constructor.

      // Compiler destructor.
  
      /**
      * Register an observer.
      *
      * \param observer  observer object (invokes method)
      * \param methodPtr pointer to relevant method
      */
      template <class Observer>
      void addObserver(Observer& observer, void (Observer::*methodPtr)(const T&));
   
      /**
      * Notify all observers.
      */
      void notify(const T& t);
   
   private:
   
      /// A linked list of functors associated with member functions.
      std::list<IFunctor<T>*> functorPtrs_;

   };

   /* 
   * Register an observer (add to list).
   */
   template <typename T>
   template <class Observer> void 
   Signal<T>::addObserver(Observer& observer, void (Observer::*methodPtr)(const T&))
   {  functorPtrs_.push_back(new MethodFunctor<Observer, T>(&observer, methodPtr)); }

   /* 
   * Notify observers (call associated methods).
   */
   template <typename T>
   void Signal<T>::notify(const T& t)
   {
      typename std::list< IFunctor<T>* >::iterator pos;
      pos = functorPtrs_.begin();
      while (pos != functorPtrs_.end())
      {
         (**pos)(t);
         ++pos;
      }
   }

   // Signal with no arguments.

   /**
   * Notifier (or subject) in the Observer design pattern.
   *
   * A Signal manages a list of registered functor objects, and provides
   * a notify() method that calls them all with the same argument.
   * 
   * \ingroup Util_Module
   */
   template <>
   class Signal<void>
   {
   
   public:

      // Compiler default constructor.

      // Compiler destructor.
  
      /**
      * Register an observer.
      *
      * \param observer  observer object (invokes method)
      * \param methodPtr pointer to relevant method
      */
      template <class Observer>
      void addObserver(Observer& observer, void (Observer::*methodPtr)());
   
      /**
      * Notify all observers.
      */
      void notify();
   
   private:
   
      /// A linked list of functors associated with member functions.
      std::list<IFunctor<>*> functorPtrs_;

   };

   /* 
   * Register an observer (add to list).
   */
   template <class Observer> void 
   Signal<>::addObserver(Observer& observer, void (Observer::*methodPtr)())
   {  functorPtrs_.push_back(new MethodFunctor<Observer>(&observer, methodPtr)); }

   /* 
   * Notify observers (call associated methods).
   */
   void Signal<>::notify()
   {
      typename std::list< IFunctor<>* >::iterator pos;
      pos = functorPtrs_.begin();
      while (pos != functorPtrs_.end())
      {
         (**pos)();
         ++pos;
      }
   }

}
#endif 
