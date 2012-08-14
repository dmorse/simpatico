#ifndef DDMD_ACTION_FACTORY_H
#define DDMD_ACTION_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>    // base class template
#include "Action.h"                // base template parameter

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Factory for DdMd::Action objects.
   *
   * \ingroup DdMd_Factory_Module
   * \ingroup DdMd_Action_Module
   */
   class ActionFactory : public Factory<Action>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation     parent Simulation
      */
      ActionFactory(Simulation& simulation);

      /** 
      * Return pointer to a new Action object.
      *
      * \param  className name of a subclass of Action.
      * \return base class pointer to a new instance of className.
      */
      virtual Action* factory(const std::string& className) const;

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
