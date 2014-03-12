#ifndef DDMD_MODIFIER_FACTORY_H
#define DDMD_MODIFIER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>    // base class template
#include "Modifier.h"                // base template parameter

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Factory for DdMd::Modifier objects.
   *
   * \ingroup DdMd_Factory_Module
   * \ingroup DdMd_Modifier_Module
   */
   class ModifierFactory : public Factory<Modifier>
   {

   public:

      /**
      * Default constructor (for unit testing)
      */
      ModifierFactory();

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object
      */
      ModifierFactory(Simulation& simulation);

      /** 
      * Return pointer to a new Modifier object.
      *
      * \param  className name of a subclass of Modifier.
      * \return base class pointer to a new instance of className.
      */
      virtual Modifier* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent Simulation.
      */
      Simulation& simulation() const
      {
         assert(simulationPtr_);  
         return *simulationPtr_; 
      }

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

   };

}
#endif
