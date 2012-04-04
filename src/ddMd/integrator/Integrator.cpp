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
     : SimulationAccess(simulation),
       timer_(NTime)
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
      timer().reduce(domain().communicator());
      atomStorage().computeNAtomTotal(domain().communicator());
      pairPotential().pairList().computeStatistics(domain().communicator());

      if (domain().isMaster()) {

         int nAtomTot = atomStorage().nAtomTotal();
         int nProc = 1;
         #ifdef UTIL_MPI
         nProc = domain().communicator().Get_size();
         #endif

         double time = timer().time();;

         // Output total time for the run
         out << std::endl;
         out << "Time Statistics" << std::endl;
         out << "nStep                " << nStep_ << std::endl;
         out << "run time             " << time << " sec" << std::endl;
         out << "time / nStep         " << time/double(nStep_) 
             << " sec" << std::endl;

         double ratio = double(nProc)/double(nStep_*nAtomTot);
         double diagnosticT =  timer().time(DIAGNOSTIC);
         double integrate1T =  timer().time(INTEGRATE1);
         double exchangeT =  timer().time(EXCHANGE);
         double neighborT =  timer().time(NEIGHBOR);
         double updateT =  timer().time(UPDATE);
         double forceT =  timer().time(FORCE);
         double integrate2T =  timer().time(INTEGRATE2);

         out << std::endl;
         out << "time * nproc / (nStep*nAtom):" << std::endl;
         out << "Total                " << Dbl(time*ratio, 12, 6)
             << " sec   " << std::endl;
         out << "Diagnostics          " << Dbl(diagnosticT*ratio, 12, 6)
             << " sec   " << Dbl(diagnosticT/time, 12, 6, true) << std::endl;
         out << "Integrate1           " << Dbl(integrate1T*ratio, 12, 6) 
             << " sec   " << Dbl(integrate1T/time, 12, 6, true) << std::endl;
         out << "Exchange             " << Dbl(exchangeT*ratio, 12, 6)
             << " sec   " << Dbl(exchangeT/time, 12, 6, true) << std::endl;
         out << "Neighbor             " << Dbl(neighborT*ratio, 12, 6)
             << " sec   " << Dbl(neighborT/time, 12, 6, true) << std::endl;
         out << "Update               " << Dbl(updateT*ratio, 12, 6)
             << " sec   " << Dbl(updateT/time, 12, 6, true) << std::endl;
         out << "Force                " << Dbl(forceT*ratio, 12, 6)
             << " sec   " << Dbl(forceT/time, 12 , 6, true) << std::endl;
         out << "Integrate2           " << Dbl(integrate2T*ratio, 12, 6) 
             << " sec   " << Dbl(integrate2T/time, 12, 6, true) << std::endl;
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
