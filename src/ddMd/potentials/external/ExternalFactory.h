#ifndef DDMD_EXTERNAL_FACTORY_H
#define DDMD_EXTERNAL_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                  // base class template
#include <ddMd/potentials/external/ExternalPotential.h>  // template argument

#include <string>

namespace DdMd
{

   class Simulation;

   /**
   * Factory for ExternalPotential objects.
   *
   * \ingroup External_Module
   */
   class ExternalFactory : public Factory<ExternalPotential>
   {

   public:
   
      /**
      * Default constructor.
      */
      ExternalFactory(Simulation& simulation);

      /**
      * Return a pointer to a new ExternalInteration, if possible.
      */
      ExternalPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to the parent Simulation.
      Simulation* simulationPtr_;

   };
  
}
#endif
