#ifndef MCMD_MC_MODULE_H
#define MCMD_MC_MODULE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <vector>

namespace McMd
{

   class McSimulation;
   class McSystem;

   /**
   * Abstract Module for an McSimulation.
   *
   * An McModule organizes a set of factory classes for an McSimulation.
   * Usage: In the main program, create an instance of the desired subclass
   * of McModule, then call addFactories().
   *
   * \ingroup McMd_Simulation_Module
   */
   class McModule
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation the associated McSimulation.
      */
      McModule(McSimulation& simulation);

      /**
      * Destructor.
      */
      virtual ~McModule();

      /**
      * Adds all factories of all types to the simulation.
      */
      virtual void addFactories() = 0;

   protected:

      /**
      * Return McSimulation by reference.
      */
      McSimulation& simulation() const;

      /**
      * Return McSystem by reference.
      */
      McSystem& system() const;

   private:

      // Pointer to parent McSimulation.
      McSimulation* simulationPtr_;

      // Pointer to parent McSystem.
      McSystem* systemPtr_;

   };

}
#endif
