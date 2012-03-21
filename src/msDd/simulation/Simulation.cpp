#ifndef MSDD_SIMULATION_CPP
#define MSDD_SIMULATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <msDd/simulation/Simulation.h>       // class header
#include <mcMd/mcSimulation/McSimulation.h>                    
#include <ddMd/system/System.h>       
//#include <ddMd/communicator/GroupCollector_inc.h>

namespace MsDd
{

   /*
   * Constructor.
   */
   Simulation::Simulation(MPI::Intracomm& ddCommunicator)
    : ddCommunicatorPtr_(&ddCommunicator),
      mcSimulationPtr_(0),
      ddSystemPtr_(0),
      isDdMaster_(false)
   {
      int myRank = ddCommunicator.Get_rank();
      if (myRank == 0) {
         isDdMaster_ = true;
         mcSimulationPtr_ = new McMd::McSimulation();
      }
   }

   /*
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
   void Simulation::readParam(std::istream& in)
   {
   }

   /*
   * Read and broadcast commands.
   *
   * Call only on master processor. Implements main
   * command loop for master processor.
   */
   void Simulation::readCommands()
   {}

   /*
   * Receive commands from master processor.
   *
   * Call only on slave processor. Allows master to control
   * action of slaves. Implements main loop for slaves.
   */
   void Simulation::receiveCommands()
   {}

   /*
   * Distribute atoms and groups from master to slaves.
   */
   void Simulation::distribute()
   {}

   /*
   * Collect atoms from master to slaves.
   */
   void Simulation::collect()
   {}


}
#endif
