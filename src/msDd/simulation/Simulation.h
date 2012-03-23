#ifndef MSDD_SIMULATION_H
#define MSDD_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>           // base class

//#include <ddMd/communicate/AtomDistributor.h>    // member
//#include <ddMd/communicate/GroupDistributor.h>   // member
//#include <ddMd/communicate/AtomCollector.h>      // member
//#include <ddMd/communicate/GroupCollector.h>     // member


namespace McMd {
   class McSimulation;
}

namespace DdMd {
   class System;
}

namespace MsDd
{

   using namespace Util;

   /**
   * Main object in master-slave simulation.
   *
   * A MsDd::Simulation object is the main object in a simulation
   * in which a McMd::McSimulation MC simulation, which exists 
   * only on the master node, controls a parallel DdMd::System
   * molecular dynamics simulation.
   *
   * Usage:
   *
   * Simulation sim(MPI::WORLD_COMM);
   * sim.readParam(std::cin);
   * if (sim.isMaster()) {
   *    sim.readCommands();
   * } else {
   *    sim.receiveCommands();
   * }
   */
   class Simulation : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Simulation(MPI::Intracomm& ddCommunicator = MPI::COMM_WORLD);

      /**
      * Destructor.
      */
      ~Simulation();

      /**
      * Initialize parent McMd::McSimulation and DdMd::System
      *
      * On master:
      *    Call McMd::McSimulation::readParam().
      *    Generate param file for DdMd::System.
      *    Send signal to read file.
      *    Call DdMd::System::readParam().
      * On all processors:
      *    Receive signal to read file.
      *    Call DdMd::System::readParam().
      */
      virtual void readParam(std::istream& in);

      /**
      * Read and broadcast commands.
      *
      * Call only on master processor. Implements main
      * command loop for master processor.
      */
      void readCommands();

      /**
      * Receive commands from master processor.
      *
      * Call only on slave processor. Allows master to control
      * action of slaves. Implements main loop for slaves.
      */
      void receiveCommands();

      /**
      * Distribute atoms and groups from master to slaves.
      */
      void distribute();

      /**
      * Collect atoms from master to slaves.
      */
      void collect();

      /**
      * Return McMd::McSimulation (valid only on master)
      */
      McMd::McSimulation& mcSimulation();

      /**
      * Return DdMd::System (valid on any processor)
      */
      DdMd::System& ddMdSystem();

      /**
      * Is this the master node in the intracommunicator?
      */
      bool isDdMaster() const;

   private:

      /// Parent McSimulation (exists only on master).
      McMd::McSimulation* mcSimulationPtr_;

      /// Slave parallel MD System (exists on all modes).
      DdMd::System*       ddSystemPtr_;

      //DdMd::AtomDistributor     atomDistributor_;
      //DdMd::GroupDistributor<2> bondDistributor_;

      //DdMd::AtomCollector       atomCollector_;
      //DdMd::GroupCollector<2>   bondCollector_;

      /// Intracommunicator for DdMd simulation.
      MPI::Intracomm* ddCommunicatorPtr_;

      /// Is this the master (rank = 0) node of the ddCommunicator?
      bool isDdMaster_;

   };

   /*
   * Return McMd::McSimulation (valid only on master)
   */
   inline
   McMd::McSimulation& Simulation::mcSimulation()
   {   return *mcSimulationPtr_; }

   /**
   * Return DdMd::System (valid on any processor)
   */
   inline 
   DdMd::System& Simulation::ddMdSystem()
   {  return *ddSystemPtr_; }

   /**
   * Is this the master node in the intracommunicator?
   */
   inline 
   bool Simulation::isDdMaster() const
   {  return isDdMaster_; }

}
#endif
