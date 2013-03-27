#ifndef DDMD_INTEGRATOR_CPP
#define DDMD_INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Integrator.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/potentials/bond/BondPotential.h>
#ifdef INTER_ANGLE
#include <ddMd/potentials/angle/AnglePotential.h>
#endif
#ifdef INTER_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef INTER_EXTERNAL
#include <ddMd/potentials/external/ExternalPotential.h>
#endif
#include <util/ensembles/BoundaryEnsemble.h>

#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/format/Bool.h>
#include <util/global.h>

namespace DdMd
{

   class Simulation;

   /*
   * Constructor.
   */
   Integrator::Integrator(Simulation& simulation)
     : SimulationAccess(simulation),
       timer_(Integrator::NTime),
       isSetup_(false)
   {}

   /*
   * Destructor.
   */
   Integrator::~Integrator()
   {}

   /*
   * Initialize atom distribution, AtomStorage, PairList and forces.
   */
   void Integrator::setupAtoms()
   {
      atomStorage().clearSnapshot();
      exchanger().exchange();
      pairPotential().buildCellList();
      if (!UTIL_ORTHOGONAL) {
         atomStorage().transformGenToCart(boundary());
      }
      atomStorage().makeSnapshot();
      pairPotential().buildPairList();
      if (simulation().boundaryEnsemble().isRigid()) {
         simulation().computeForces();
      } else {
         simulation().computeForcesAndVirial();
      }
   }

   /*
   * Compute forces for all atoms, with timing.
   */
   void Integrator::computeForces()
   {
      simulation().zeroForces();
      pairPotential().computeForces();
      timer_.stamp(PAIR_FORCE);
      bondPotential().computeForces();
      timer_.stamp(BOND_FORCE);
      #ifdef INTER_ANGLE
      if (nAngleType()) {
         anglePotential().computeForces();
         timer_.stamp(ANGLE_FORCE);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType()) {
         dihedralPotential().computeForces();
         timer_.stamp(DIHEDRAL_FORCE);
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal()) {
         externalPotential().computeForces();
      }
      #endif

      // Reverse communication (if any)
      if (reverseUpdateFlag()) {
         exchanger().reverseUpdate();
      }

      // Send signal indicating change in atomic forces
      // simulation().forceSignal().notify();
   }

   /*
   * Compute forces for all local atoms and virial, with timing.
   */
   void Integrator::computeForcesAndVirial()
   {
      simulation().zeroForces();
      pairPotential().computeForcesAndStress(domain().communicator());
      timer_.stamp(PAIR_FORCE);
      bondPotential().computeForcesAndStress(domain().communicator());
      timer_.stamp(BOND_FORCE);
      #ifdef INTER_ANGLE
      if (nAngleType()) {
         anglePotential().computeForcesAndStress(domain().communicator());
         timer_.stamp(ANGLE_FORCE);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType()) {
         dihedralPotential().computeForcesAndStress(domain().communicator());
         timer_.stamp(DIHEDRAL_FORCE);
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal()) {
         externalPotential().computeForces();
      }
      #endif

      // Reverse communication (if any)
      if (reverseUpdateFlag()) {
         exchanger().reverseUpdate();
      }

      // Send signal indicating change in atomic forces
      // simulation().forceSignal().notify();
   }

   /*
   * Determine whether an atom exchange and reneighboring is needed.
   */
   bool Integrator::isExchangeNeeded(double skin) 
   {
     
      if (!atomStorage().isCartesian()) {
         UTIL_THROW("Error: Coordinates not Cartesian in isExchangeNeeded");
      } 

      // Calculate maximum square displacment on this node
      double maxSqDisp = atomStorage().maxSqDisplacement(); 
      int    needed = 0;
      if (sqrt(maxSqDisp) > 0.5*skin) {
         needed = 1; 
      }
      timer_.stamp(CHECK);

      #if UTIL_MPI
      int neededAll;
      domain().communicator().Allreduce(&needed, &neededAll, 
                                        1, MPI::INT, MPI::MAX);
      timer_.stamp(ALLREDUCE);
      return bool(neededAll);
      #else
      return bool(needed);
      #endif
   }

   #if 0
   /*
   * Determine whether an atom exchange and reneighboring is needed.
   */
   bool Integrator::isExchangeNeeded(double skin) 
   {
     
      if (!atomStorage().isCartesian()) {
         UTIL_THROW("Error: Coordinates not Cartesian in isExchangeNeeded");
      } 

      // Calculate maximum square displacment on this node
      double maxSqDisp = atomStorage().maxSqDisplacement(); 
      timer_.stamp(CHECK);

      // Decide on master node if maximum exceeds threshhold.
      int needed;

      #if UTIL_MPI
      double maxSqDispAll;                    // global maximum
      domain().communicator().Reduce(&maxSqDisp, &maxSqDispAll, 1, 
                          MPI::DOUBLE, MPI::MAX, 0);
      if (domain().communicator().Get_rank() == 0) {
         needed = 0;
         if (sqrt(maxSqDispAll) > 0.5*skin) {
            needed = 1; 
         }
      }
      domain().communicator().Bcast(&needed, 1, MPI::INT, 0);
      timer_.stamp(ALLREDUCE);
      #else
      if (sqrt(maxSqDisp) > 0.5*skin) {
         needed = 1; 
      }
      #endif

      return bool(needed);
   }
   #endif

   /*
   * Return time per processor for last run.
   */
   double Integrator::time() const
   {  return timer_.time(); }

   /*
   * Reduce timing statistics data from all processors.
   */
   void Integrator::computeStatistics()
   {  
      #ifdef UTIL_MPI
      timer().reduce(domain().communicator());  
      #endif
   }

   /*
   * Output statistics.
   */
   void Integrator::outputStatistics(std::ostream& out)
   {
      if (!domain().isMaster()) {
         UTIL_THROW("May be called only on domain master");
      }

      double time = timer().time();
      int nAtomTot = atomStorage().nAtomTotal();
      int nProc = 1;
      #ifdef UTIL_MPI
      nProc = domain().communicator().Get_size();
      #endif

      // Output total time for the run
      out << std::endl;
      out << "Time Statistics" << std::endl;
      out << "nStep                " << iStep_ << std::endl;
      out << "run time             " << time << " sec" << std::endl;
      out << "time / nStep         " << time/double(iStep_) 
          << " sec" << std::endl;


      double factor1 = 1.0/double(iStep_);
      double factor2 = double(nProc)/double(iStep_*nAtomTot);
      double totalT = 0.0;

      out << std::endl;
      out << "T = Time per processor, M = nstep = # steps" << std::endl
          << "P = # procs, N = # atoms (total, all processors)"
          << std::endl << std::endl;
      out << "                     " 
          << "   T/M [sec]   "
          << "   T*P/(N*M)   "
          << " Percent (%)" << std::endl;
      out << "Total                " 
          << Dbl(time*factor1, 12, 6)
          << "   "
          << Dbl(time*factor2, 12, 6)
          << "   " << Dbl(100.0, 12, 6, true) << std::endl;
      double diagnosticT = timer().time(DIAGNOSTIC);
      totalT += diagnosticT;
      out << "Diagnostics          " 
          << Dbl(diagnosticT*factor1, 12, 6)
          << "   "
          << Dbl(diagnosticT*factor2, 12, 6)
          << "   " << Dbl(100.0*diagnosticT/time, 12, 6, true) << std::endl;
      double integrate1T = timer().time(INTEGRATE1);
      totalT += integrate1T;
      out << "Integrate1           " 
          << Dbl(integrate1T*factor1, 12, 6) 
          << "   "
          << Dbl(integrate1T*factor2, 12, 6) 
          << "   " << Dbl(100.0*integrate1T/time, 12, 6, true) << std::endl;
      double checkT =  timer().time(CHECK);
      totalT += checkT;
      out << "Check                " 
          << Dbl(checkT*factor1, 12, 6)
          << "   "
          << Dbl(checkT*factor2, 12, 6)
          << "   " << Dbl(100.0*checkT/time, 12, 6, true) << std::endl;
      double allReduceT =  timer().time(ALLREDUCE);
      totalT += allReduceT;
      out << "AllReduce            " 
          << Dbl(allReduceT*factor1, 12, 6)
          << "   "
          << Dbl(allReduceT*factor2, 12, 6)
          << "   " << Dbl(100.0*allReduceT/time, 12, 6, true) << std::endl;
      if (!UTIL_ORTHOGONAL) {
         double transformFT = timer().time(TRANSFORM_F);
         totalT += transformFT;
         out << "Transform (forward)  " 
             << Dbl(transformFT*factor1, 12, 6)
          << "   "
             << Dbl(transformFT*factor2, 12, 6)
             << "   " << Dbl(100.0*transformFT/time, 12, 6, true) << std::endl;
      }
      double exchangeT = timer().time(EXCHANGE);
      totalT += exchangeT;
      out << "Exchange             " 
          << Dbl(exchangeT*factor1, 12, 6)
          << "   "
          << Dbl(exchangeT*factor2, 12, 6)
          << "   " << Dbl(100.0*exchangeT/time, 12, 6, true) << std::endl;
      double cellListT = timer().time(CELLLIST);
      totalT += cellListT;
      out << "CellList             " 
          << Dbl(cellListT*factor1, 12, 6)
          << "   "
          << Dbl(cellListT*factor2, 12, 6)
          << "   " << Dbl(100.0*cellListT/time, 12, 6, true) << std::endl;
      if (!UTIL_ORTHOGONAL) {
         double transformRT = timer().time(TRANSFORM_R);
         totalT += transformRT;
         out << "Transform (reverse)  " 
             << Dbl(transformRT*factor1, 12, 6)
             << "   "
             << Dbl(transformRT*factor2, 12, 6)
             << "   " << Dbl(100.0*transformRT/time, 12, 6, true) << std::endl;
      }
      double pairListT = timer().time(PAIRLIST);
      totalT += pairListT;
      out << "PairList             " 
          << Dbl(pairListT*factor1, 12, 6)
          << "   "
          << Dbl(pairListT*factor2, 12, 6)
          << "   " << Dbl(100.0*pairListT/time, 12, 6, true) << std::endl;
      double updateT = timer().time(UPDATE);
      totalT += updateT;
      out << "Update               " 
          << Dbl(updateT*factor1, 12, 6)
          << "   "
          << Dbl(updateT*factor2, 12, 6)
          << "   " << Dbl(100.0*updateT/time, 12, 6, true) << std::endl;
      double pairForceT = timer().time(PAIR_FORCE);
      totalT += pairForceT;
      out << "Pair Forces          " 
          << Dbl(pairForceT*factor1, 12, 6)
          << "   "
          << Dbl(pairForceT*factor2, 12, 6)
          << "   " << Dbl(100.0*pairForceT/time, 12 , 6, true) << std::endl;
      double bondForceT = timer().time(BOND_FORCE);
      totalT += bondForceT;
      out << "Bond Forces          " 
          << Dbl(bondForceT, 12, 6)
          << "   "
          << Dbl(bondForceT*factor2, 12, 6)
          << "   " << Dbl(100.0*bondForceT/time, 12 , 6, true) << std::endl;
      #ifdef INTER_ANGLE
      double angleForceT = timer().time(ANGLE_FORCE);
      totalT += angleForceT;
      out << "Angle Forces         " 
          << Dbl(angleForceT, 12, 6)
          << "   "
          << Dbl(angleForceT*factor2, 12, 6)
          << "   " << Dbl(100.0*angleForceT/time, 12 , 6, true) << std::endl;
      #endif
      #ifdef INTER_DIHEDRAL
      double dihedralForceT = timer().time(DIHEDRAL_FORCE);
      totalT += dihedralForceT;
      out << "Dihedral Forces      " 
          << Dbl(dihedralForceT, 12, 6) 
          << "   "
          << Dbl(dihedralForceT*factor2, 12, 6)
          << "   " << Dbl(100.0*dihedralForceT/time, 12 , 6, true) << std::endl;
      #endif
      double integrate2T = timer().time(INTEGRATE2);
      totalT += integrate2T;
      out << "Integrate2           " 
          << Dbl(integrate2T, 12, 6) 
          << "   "
          << Dbl(integrate2T*factor2, 12, 6) 
          << "   " << Dbl(100.0*integrate2T/time, 12, 6, true) << std::endl;
      //out << "Sum of above         " << Dbl(totalT*ratio, 12, 6) 
      //    << " sec   " << Dbl(100.0*totalT/time, 12, 6, true) << std::endl;
      out << std::endl;

      int buildCounter = pairPotential().pairList().buildCounter(); 
      out << "buildCounter             " 
                  << Int(buildCounter, 10)
                  << std::endl;
      out << "steps / build            "
                  << double(iStep_)/double(buildCounter)
                  << std::endl;
      out << std::endl;

   }

   /**
   * Clear timing, dynamical state, statistics, and diagnostic accumulators.
   */
   void Integrator::clear()
   { 
      iStep_ = 0;
      initDynamicalState();
      timer().clear(); 
      simulation().exchanger().timer().clear();
      simulation().buffer().clearStatistics();
      pairPotential().pairList().clearStatistics();
      atomStorage().clearStatistics();
      bondStorage().clearStatistics();
      #ifdef INTER_ANGLE
      angleStorage().clearStatistics();
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage().clearStatistics();
      #endif
      simulation().diagnosticManager().clear();
   }

}
#endif
