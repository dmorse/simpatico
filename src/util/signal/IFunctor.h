#ifndef I_FUNCTOR_H
#define I_FUNCTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <list>

namespace Util
{

   /**
   * Interface for functor that wraps a void function with one argument (abstract).
   */
   template <typename T=void>
   class IFunctor {

   public:

      /**
      * Destructor (virtual)
      */
      virtual ~IFunctor(){}

      /**
      * Call the associated function.
      *
      * \param t parameter value passed to associated function.
      */
      virtual void operator () (const T&) = 0;

   };

   // Functor for functions with no arguments.

   /**
   * Interface for functor that wraps a void function with no arguments (abstract).
   *
   * The operator (const T&) invokes the associated zero-parameter function.
   */
   template <>
   class IFunctor<void> {

   public: 

      /**
      * Destructor (virtual)
      */
      virtual ~IFunctor(){}

      /**
      * Call a specific member function with one parameter.
      */
      virtual void operator () () = 0; 

   };

}
#endif 
