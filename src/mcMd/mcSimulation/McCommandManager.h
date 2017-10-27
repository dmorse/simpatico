#ifndef MCMD_MC_COMMAND_MANAGER_H
#define MCMD_MC_COMMAND_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/commands/CommandManager.h> // base class 

namespace McMd
{

   class McSimulation;
   class McSystem;

   using namespace Util;

   /**
   * Command interpreter and Manager for an McSimulation.
   *
   * The implementation of the readStandardCommand() function defines the 
   * standard, built-in commands for an McSimulation. Additional commands
   * may be added in the optional CommandManager section of the parameter
   * file and then invoked in the command file. 
   *
   * The default Factory<Command> object for an McCommandManager is an
   * McCommandFactory.
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

      /**
      * Get parent McSimulation by reference.
      */
      McSimulation& simulation() const;

      /**
      * Get associated McSystem by reference.
      */
      McSystem& system() const;

      /**
      * Attempt to read one of the standard commands.
      *
      * Precondition: The command identifier should have been read from
      * stream in, and should be passed as string parameter command. 
      * If the command name is matched, any further parameters required 
      * by the command are read from stream in.
      *
      * \param command command name (capitalized identifer)
      * \param in input stream from which command was read
      */
      virtual
      bool readStandardCommand(std::string command, std::istream& in);

   private:

      /// Pointer to parent Simulation
      McSimulation*  simulationPtr_;

      /// Pointer to associated McSystem.
      McSystem*  systemPtr_;

   };


   // Inline functions

   inline
   McSimulation& McCommandManager::simulation() const
   {
      assert(simulationPtr_);
      return *simulationPtr_;
   }

   inline
   McSystem& McCommandManager::system() const
   {
      assert(systemPtr_);
      return *systemPtr_;
   }

}
#endif
