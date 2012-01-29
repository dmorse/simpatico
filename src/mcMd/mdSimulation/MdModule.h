#ifndef MD_MODULE_H
#define MD_MODULE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class MdSimulation;
   class MdSystem;

   /**
   * A Module for an MdSimulation.
   *
   * \ingroup Simulation_Module
   */
   class MdModule
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent MdSimulation
      */
      MdModule(MdSimulation& simulation);

      /**
      * Destructor.
      */
      virtual ~MdModule();

      /**
      * Adds all factories of all types to the simulation.
      */
      virtual void addFactories() = 0;

   protected:

      /**
      * Return MdSimulation by reference.
      */
      MdSimulation& simulation() const;

      /**
      * Return MdSystem by reference.
      */
      MdSystem& system() const;

   private:

      // Pointer to parent MdSimulation.
      MdSimulation* simulationPtr_;

      // Pointer to parent MdSystem.
      MdSystem* systemPtr_;

   };

}
#endif
