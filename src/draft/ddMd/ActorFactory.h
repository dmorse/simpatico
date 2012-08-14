#ifndef DDMD_ACTOR_FACTORY_H
#define DDMD_ACTOR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>    // base class template
#include "Actor.h"                // base template parameter

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Factory for DdMd::Actor objects.
   *
   * \ingroup DdMd_Factory_Module
   * \ingroup DdMd_Actor_Module
   */
   class ActorFactory : public Factory<Actor>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation     parent Simulation
      */
      ActorFactory(Simulation& simulation);

      /** 
      * Return pointer to a new Actor object.
      *
      * \param  className name of a subclass of Actor.
      * \return base class pointer to a new instance of className.
      */
      virtual Actor* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent Simulation.
      */
      Simulation& simulation() const
      {  return *simulationPtr_; }

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

   };

}
#endif
