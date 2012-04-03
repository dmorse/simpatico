#ifndef DDMD_INTEGRATOR_CPP
#define DDMD_INTEGRATOR_CPP

#include "Integrator.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/potentials/pair/PairPotential.h>

#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/format/Bool.h>
#include <util/util/Timer.h>
#include <util/global.h>


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class Simulation;

   /*
   * Constructor.
   */
   Integrator::Integrator(Simulation& simulation)
     : SimulationAccess(simulation)
   {}

   /*
   * Destructor.
   */
   Integrator::~Integrator()
   {}

   /*
   * Output statistics.
   */
   void Integrator::outputStatistics(std::ostream& out)
   {
      atomStorage().computeNAtomTotal(domain().communicator());
      pairPotential().pairList().computeStatistics(domain().communicator());

      if (domain().isMaster()) {

         int nAtomTot = atomStorage().nAtomTotal();
         int nProc = 1;
         #ifdef UTIL_MPI
         nProc = domain().communicator().Get_size();
         #endif

         // Output total time for the run
         out << std::endl;
         out << "Time Statistics" << std::endl;
         out << "nStep                " << nStep_ << std::endl;
         out << "run time             " << timer().time() << " sec" << std::endl;
         out << "time / nStep         " << timer().time()/double(nStep_) 
   	             << " sec" << std::endl;
         out << "time / (nStep*nAtom) " 
                     << timer().time()*double(nProc)/double(nStep_*nAtomTot)
                     << " sec" << std::endl;
         out << std::endl;
         out << std::endl;
 
         //pairPotential().pairList().outputStatistics(out);
         out << "PairList Statistics" << std::endl;
         out << "maxNPair, capacity " 
                     << Int(pairPotential().pairList().maxNPair(), 10)
                     << Int(pairPotential().pairList().pairCapacity(), 10)
                     << std::endl;
         out << "maxNAtom, capacity " 
                     << Int(pairPotential().pairList().maxNAtom(), 10)
                     << Int(pairPotential().pairList().atomCapacity(), 10)
                     << std::endl;
         out << "buildCounter       " 
                     << Int(pairPotential().pairList().buildCounter(), 10)
                     << std::endl;
         out << "steps / build      "
                     << double(nStep_)/double(pairPotential().pairList().buildCounter())
                     << std::endl;
         out << std::endl;

      }

   }

}
#endif
