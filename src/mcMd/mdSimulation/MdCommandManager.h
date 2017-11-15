#ifndef MCMD_MD_COMMAND_MANAGER_H
#define MCMD_MD_COMMAND_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/commands/CommandManager.h>   // base class

namespace Util { template <typename T> class Factory; }

namespace McMd
{

   class Simulation;
   class MdSimulation;
   class MdSystem;

   using namespace Util;

   /**
   * Command interpreter and Command Manager for an MdSimulation.
   *
   * The implementation of the readStandardCommand() function defines the 
   * standard, built-in commands for an MdSimulation. Additional commands
   * may be added in the optional CommandManager section of the parameter
   * file and then invoked in the command file. 
   *
   * The default Factory<Command> object for an MdCommandManager is an
   * MdCommandFactory.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_Command_Md_Module
   */
   class MdCommandManager : public CommandManager
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent MdSimulation
      */
      MdCommandManager(MdSimulation& simulation);

      /**
      * Constructor.
      *
      * \param simulation parent Simulation
      * \param system     associated MdSystem
      */
      MdCommandManager(MdSimulation& simulation, MdSystem& system);

      /**
      * Destructor.
      */
      virtual ~MdCommandManager();

   protected:

      /**
      * Create and return pointer to a new MdCommandFactory object
      */
      virtual Util::Factory<Command>* newDefaultFactory() const;

      /**
      * Get parent MdSimulation by reference.
      */
      MdSimulation& simulation() const;

      /**
      * Get associated MdSystem by reference.
      */
      MdSystem& system() const;

      /**
      * Attempt to read one of the standard built-in commands.
      *
      * Precondition: The command identifier should have been read from
      * stream in, and should be passed as string parameter command. 
      * If the command name is matched, any further parameters required 
      * by the command are read from stream in, the command is executed,
      * the function returns a true value. If the command name is not
      * recognized, this function immediately returns false.
      *
      * \param command command name (capitalized identifer)
      * \param in input stream from which command was read
      */
      virtual
      bool readStandardCommand(std::string command, std::istream& in);

   private:

      /// Pointer to parent Simulation
      MdSimulation* simulationPtr_;

      /// Pointer to associated MdSystem.
      MdSystem*   systemPtr_;

   };

   // Inline member functions

   inline
   MdSimulation& MdCommandManager::simulation() const
   {
      assert(simulationPtr_);
      return *simulationPtr_;
   }

   inline
   MdSystem& MdCommandManager::system() const
   {
      assert(systemPtr_);
      return *systemPtr_;
   }

}
#endif
