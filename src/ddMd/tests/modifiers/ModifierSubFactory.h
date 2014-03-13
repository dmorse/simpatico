#ifndef DDMD_MODIFIER_SUB_FACTORY_H
#define DDMD_MODIFIER_SUB_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/modifiers/ModifierFactory.h>  // base class 
#include "ModifierClasses.h"                 // Subclasses for factory

namespace DdMd
{

   using namespace Util;

   // class Simulation;

   /*
   * Factory for DdMd::Modifier objects.
   */
   class ModifierSubFactory : public ModifierFactory
   {

   public:

      /*
      * Default constructor (for unit testing).
      */
      ModifierSubFactory();

      /*
      * Constructor.
      *
      * \param simulation     parent Simulation
      */
      ModifierSubFactory(Simulation& simulation);

      /*
      * Return pointer to a new Modifier object.
      *
      * \param  className name of a subclass of Modifier.
      * \return base class pointer to a new instance of className.
      */
      virtual Modifier* factory(const std::string& className) const;

   };

   /*
   * Constructor.
   */
   ModifierSubFactory::ModifierSubFactory()
    : ModifierFactory()
   {}

   /* 
   * Return a pointer to an instance of Modifier subclass className.
   */
   Modifier* ModifierSubFactory::factory(const std::string &className) 
   const
   {
      Modifier* ptr = 0;

      // Try subfactories first (if any)
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Simulation Modifiers
      if (className == "ModifierA") {
         ptr = new ModifierA();
      } // else 

      return ptr;
   }

}
#endif
