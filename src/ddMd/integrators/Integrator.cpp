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
       timer_(Integrator::NTime)
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
      #ifdef UTIL_MPI
      timer().reduce(domain().communicator());
      exchanger().timer().reduce(domain().communicator());
      atomStorage().computeNAtomTotal(domain().communicator());
      pairPotential().pairList().computeStatistics(domain().communicator());
      #else
      pairPotential().pairList().computeStatistics();
      #endif

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
         double checkT =  timer().time(CHECK);
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
         out << "Check                " << Dbl(checkT*ratio, 12, 6)
             << " sec   " << Dbl(checkT/time, 12, 6, true) << std::endl;
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

         #ifdef DDMD_EXCHANGER_TIMER
         double AtomPlanT =  exchanger().timer().time(Exchanger::ATOM_PLAN);
         out << "AtomPlan              " << Dbl(AtomPlanT*ratio, 12, 6)
             << " sec   " << Dbl(AtomPlanT/time, 12, 6, true) << std::endl;
         double InitGroupPlanT =  exchanger().timer().time(Exchanger::INIT_GROUP_PLAN);
         out << "InitGroupPlan         " << Dbl(InitGroupPlanT*ratio, 12, 6)
             << " sec   " << Dbl(InitGroupPlanT/time, 12, 6, true) << std::endl;
         double ClearGhostsT =  exchanger().timer().time(Exchanger::CLEAR_GHOSTS);
         out << "ClearGhosts           " << Dbl(ClearGhostsT*ratio, 12, 6)
             << " sec   " << Dbl(ClearGhostsT/time, 12, 6, true) << std::endl;
         double PackAtomsT =  exchanger().timer().time(Exchanger::PACK_ATOMS);
         out << "PackAtoms             " << Dbl(PackAtomsT*ratio, 12, 6)
             << " sec   " << Dbl(PackAtomsT/time, 12, 6, true) << std::endl;
         double PackGroupsT =  exchanger().timer().time(Exchanger::PACK_GROUPS);
         out << "PackGroups            " << Dbl(PackGroupsT*ratio, 12, 6)
             << " sec   " << Dbl(PackGroupsT/time, 12, 6, true) << std::endl;
         double RemoveAtomsT =  exchanger().timer().time(Exchanger::REMOVE_ATOMS);
         out << "RemoveAtoms           " << Dbl(RemoveAtomsT*ratio, 12, 6)
             << " sec   " << Dbl(RemoveAtomsT/time, 12, 6, true) << std::endl;
         double RemoveGroupsT =  exchanger().timer().time(Exchanger::REMOVE_GROUPS);
         out << "RemoveGroups          " << Dbl(RemoveGroupsT*ratio, 12, 6)
             << " sec   " << Dbl(RemoveGroupsT/time, 12, 6, true) << std::endl;
         double SendRecvAtomsT =  exchanger().timer().time(Exchanger::SEND_RECV_ATOMS);
         out << "SendRecvAtoms         " << Dbl(SendRecvAtomsT*ratio, 12, 6)
             << " sec   " << Dbl(SendRecvAtomsT/time, 12, 6, true) << std::endl;
         double UnpackAtomsT =  exchanger().timer().time(Exchanger::UNPACK_ATOMS);
         out << "UnpackAtoms           " << Dbl(UnpackAtomsT*ratio, 12, 6)
             << " sec   " << Dbl(UnpackAtomsT/time, 12, 6, true) << std::endl;
         double UnpackGroupsT =  exchanger().timer().time(Exchanger::UNPACK_GROUPS);
         out << "UnpackGroups          " << Dbl(UnpackGroupsT*ratio, 12, 6)
             << " sec   " << Dbl(UnpackGroupsT/time, 12, 6, true) << std::endl;
         double FinishGroupPlanT =  exchanger().timer().time(Exchanger::FINISH_GROUP_PLAN);
         out << "FinishGroupPlan       " << Dbl(FinishGroupPlanT*ratio, 12, 6)
             << " sec   " << Dbl(FinishGroupPlanT/time, 12, 6, true) << std::endl;
         double SendArraysT =  exchanger().timer().time(Exchanger::INIT_SEND_ARRAYS);
         out << "SendArrays            " << Dbl(SendArraysT*ratio, 12, 6)
             << " sec   " << Dbl(SendArraysT/time, 12, 6, true) << std::endl;
         double PackGhostsT =  exchanger().timer().time(Exchanger::PACK_GHOSTS);
         out << "PackGhosts            " << Dbl(PackGhostsT*ratio, 12, 6)
             << " sec   " << Dbl(PackGhostsT/time, 12, 6, true) << std::endl;
         double SendRecvGhostsT =  exchanger().timer().time(Exchanger::SEND_RECV_GHOSTS);
         out << "SendRecvGhosts        " << Dbl(SendRecvGhostsT*ratio, 12, 6)
             << " sec   " << Dbl(SendRecvGhostsT/time, 12, 6, true) << std::endl;
         double UnpackGhostsT =  exchanger().timer().time(Exchanger::UNPACK_GHOSTS);
         out << "UnpackGhosts          " << Dbl(UnpackGhostsT*ratio, 12, 6)
             << " sec   " << Dbl(UnpackGhostsT/time, 12, 6, true) << std::endl;
         double FindGroupGhostsT =  exchanger().timer().time(Exchanger::FIND_GROUP_GHOSTS);
         out << "FindGroupGhosts       " << Dbl(FindGroupGhostsT*ratio, 12, 6)
             << " sec   " << Dbl(FindGroupGhostsT/time, 12, 6, true) << std::endl;
         double PackUpdateT =  exchanger().timer().time(Exchanger::PACK_UPDATE);
         out << "PackUpdate            " << Dbl(PackUpdateT*ratio, 12, 6)
             << " sec   " << Dbl(PackUpdateT/time, 12, 6, true) << std::endl;
         double SendRecvUpdateT =  exchanger().timer().time(Exchanger::SEND_RECV_UPDATE);
         out << "SendRecvUpdate        " << Dbl(SendRecvUpdateT*ratio, 12, 6)
             << " sec   " << Dbl(SendRecvUpdateT/time, 12, 6, true) << std::endl;
         double UnpackUpdateT =  exchanger().timer().time(Exchanger::UNPACK_UPDATE);
         out << "UnpackUpdate          " << Dbl(UnpackUpdateT*ratio, 12, 6)
             << " sec   " << Dbl(UnpackUpdateT/time, 12, 6, true) << std::endl;
         double LocalUpdateT =  exchanger().timer().time(Exchanger::LOCAL_UPDATE);
         out << "LocalUpdate           " << Dbl(LocalUpdateT*ratio, 12, 6)
             << " sec   " << Dbl(LocalUpdateT/time, 12, 6, true) << std::endl;
         out << std::endl;
         out << std::endl;
         #endif

         #ifdef DDMD_PAIR_POTENTIAL_TIMER
         double BuildCellListT = pairPotential().timer().time(PairPotential::BUILD_CELL_LIST);
         out << "BuildCellList         " << Dbl(BuildCellListT*ratio, 12, 6)
             << " sec   " << Dbl(BuildCellListT/time, 12, 6, true) << std::endl;
         double BuildPairListT = pairPotential().timer().time(PairPotential::BUILD_PAIR_LIST);
         out << "BuildPairList         " << Dbl(BuildPairListT*ratio, 12, 6)
             << " sec   " << Dbl(BuildPairListT/time, 12, 6, true) << std::endl;
         double PairForcesT = pairPotential().timer().time(PairPotential::FORCES);
         out << "PairForces            " << Dbl(PairForcesT*ratio, 12, 6)
             << " sec   " << Dbl(PairForcesT/time, 12, 6, true) << std::endl;
         out << std::endl;
         out << std::endl;
         #endif

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
