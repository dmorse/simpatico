#ifndef MCMD_MC_COMMAND_FACTORY_H
#define MCMD_MC_COMMAND_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>        // base class template
#include <mcMd/commands/Command.h>     // class template param.

namespace McMd
{

   using namespace Util;

   class McSimulation;
   class McSystem;

   /**
   * CommandFactory for an McSimulation
   *
   * \ingroup McMd_Factory_Module
   * \ingroup McMd_Command_Mc_Module
   */
   class McCommandFactory : public Factory<Command>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      * \param system     parent McSystem
      */
      McCommandFactory(McSimulation& simulation, McSystem& system);

      /** 
      * Return pointer to a new Command object.
      *
      * \param  className name of a subclass of Command.
      * \return base class pointer to a new instance of className.
      */
      virtual Command* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent McSystem.
      */
      McSystem& system() const
      { return *systemPtr_; }

      /**
      * Return reference to parent McSimulation.
      */
      McSimulation& simulation() const
      { return *simulationPtr_; }

   private:

      // Pointer to parent Simulation
      McSimulation *simulationPtr_;

      // Pointer to parent System
      McSystem  *systemPtr_;

   };

}
#endif
