#ifndef MCMD_MD_COMMAND_FACTORY_H
#define MCMD_MD_COMMAND_FACTORY_H

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

   class MdSimulation;
   class MdSystem;

   /**
   * CommandFactory for an MdSimulation
   *
   * \ingroup McMd_Factory_Module
   * \ingroup McMd_Command_Module
   */
   class MdCommandFactory : public Factory<Command>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent MdSimulation
      * \param system     parent MdSystem
      */
      MdCommandFactory(MdSimulation& simulation, MdSystem& system);

      /** 
      * Return pointer to a new Command object.
      *
      * \param  className name of a subclass of Command.
      * \return base class pointer to a new instance of className.
      */
      virtual Command* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent MdSystem.
      */
      MdSystem& system() const
      { return *systemPtr_; }

      /**
      * Return reference to parent MdSimulation.
      */
      MdSimulation& simulation() const
      { return *simulationPtr_; }

   private:

      // Pointer to parent Simulation
      MdSimulation *simulationPtr_;

      // Pointer to parent System
      MdSystem  *systemPtr_;

   };

}
#endif
