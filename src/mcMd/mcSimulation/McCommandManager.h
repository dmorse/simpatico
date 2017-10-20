#ifndef MCMD_MC_COMMAND_MANAGER_H
#define MCMD_MC_COMMAND_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/commands/CommandManager.h> // base class 

namespace Util { template <typename T> class Factory; }

namespace McMd
{

   class Simulation;
   class McSimulation;
   class McSystem;

   using namespace Util;

   /**
   * Manager for Command objects in an McSimulation.
   *
   * The default Factory<Command> object for an McCommandManager is
   * an McCommandFactory.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_Command_Module
   */
   class McCommandManager : public CommandManager
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      */
      McCommandManager(McSimulation& simulation);

      /**
      * Constructor.
      *
      * \param simulation parent Simulation
      * \param system     associated McSystem
      */
      McCommandManager(McSimulation& simulation, McSystem& system);

      /**
      * Destructor.
      */
      virtual ~McCommandManager();

   protected:

      /**
      * Create and return pointer to a new McCommandFactory object
      */
      virtual Util::Factory<Command>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation
      McSimulation*  simulationPtr_;

      /// Pointer to associated McSystem.
      McSystem*  systemPtr_;

   };

}
#endif
