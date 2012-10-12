#ifndef DDMD_INTEGRATOR_CPP
#define DDMD_INTEGRATOR_CPP

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


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
       timer_(Integrator::NTime),
       isSetup_(false)
   {}

   /*
   * Destructor.
   */
   Integrator::~Integrator()
   {}

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
      simulation().forceSignal().notify();
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
      simulation().forceSignal().notify();
   }

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

      // Save contributions to total time
      double diagnosticT = timer().time(DIAGNOSTIC);
      double integrate1T = timer().time(INTEGRATE1);
      double checkT =  timer().time(CHECK);
      double transformFT = timer().time(TRANSFORM_F);
      double exchangeT = timer().time(EXCHANGE);
      double cellListT = timer().time(CELLLIST);
      double transformRT = timer().time(TRANSFORM_R);
      double pairListT = timer().time(PAIRLIST);
      double updateT = timer().time(UPDATE);
      double pairForceT = timer().time(PAIR_FORCE);
      double bondForceT = timer().time(BOND_FORCE);
      double integrate2T = timer().time(INTEGRATE2);

      // Conversion factor from time to (time per step)/(atoms per processor)
      double iStepInv = 1.0/double(iStep_); 
      double ratio    = double(nProc)/(double(iStep_)*double(nAtomTot));

      out << std::endl;
      out << "                     "
          << " T = time/nStep      " 
          << " T*nProc/nAtom       " 
          << " percentage " 
          << std::endl;
      out << "                     "
          << " ---------------     " 
          << " ---------------     " 
          << " ---------- " 
          << std::endl;
      out << "Total                " 
          << Dbl(time*iStepInv, 12, 6) << " sec     " 
          << Dbl(time*ratio, 12, 6)    << " sec   " 
          << "      .... " 
          << std::endl;
      out << "Diagnostics          " 
          << Dbl(diagnosticT*iStepInv, 12, 6) << " sec     " 
          << Dbl(diagnosticT*ratio, 12, 6)    << " sec   "
          << Dbl(diagnosticT*100.0/time, 12, 6, true) << std::endl;
      out << "Integrate1           " 
          << Dbl(integrate1T*iStepInv, 12, 6) << " sec     " 
          << Dbl(integrate1T*ratio, 12, 6)    << " sec   " 
          << Dbl(integrate1T*100.0/time, 12, 6, true) << std::endl;
      out << "Check                " 
          << Dbl(checkT*iStepInv, 12, 6) << " sec     " 
          << Dbl(checkT*ratio, 12, 6)    << " sec   " 
          << Dbl(checkT*100.0/time, 12, 6, true) << std::endl;
      if (!UTIL_ORTHOGONAL) {
         out << "Transform (forward)  " 
             << Dbl(transformFT*iStepInv, 12, 6) << " sec     " 
             << Dbl(transformFT*ratio, 12, 6)    << " sec   " 
             << Dbl(transformFT*100.0/time, 12, 6, true) << std::endl;
      }
      out << "Exchange             " 
          << Dbl(exchangeT*iStepInv, 12, 6) << " sec     " 
          << Dbl(exchangeT*ratio, 12, 6)    << " sec   " 
          << Dbl(exchangeT*100.0/time, 12, 6, true) << std::endl;
      out << "CellList             " 
          << Dbl(cellListT*iStepInv, 12, 6) << " sec     " 
          << Dbl(cellListT*ratio, 12, 6)    << " sec   " 
          << Dbl(cellListT*100.0/time, 12, 6, true) << std::endl;
      if (!UTIL_ORTHOGONAL) {
         out << "Transform (reverse)  " 
             << Dbl(transformRT*iStepInv, 12, 6) << " sec     " 
             << Dbl(transformRT*ratio, 12, 6)    << " sec   " 
             << Dbl(transformRT*100.0/time, 12, 6, true) << std::endl;
      }
      out << "PairList             " 
          << Dbl(pairListT*iStepInv, 12, 6) << " sec     " 
          << Dbl(pairListT*ratio, 12, 6)    << " sec   " 
          << Dbl(pairListT*100.0/time, 12, 6, true) << std::endl;
      out << "Update               " 
          << Dbl(updateT*iStepInv, 12, 6) << " sec     " 
          << Dbl(updateT*ratio, 12, 6)    << " sec   " 
          << Dbl(updateT*100.0/time, 12, 6, true) << std::endl;
      out << "Pair Forces          " 
          << Dbl(pairForceT*iStepInv, 12, 6) << " sec     " 
          << Dbl(pairForceT*ratio, 12, 6)    << " sec   " 
          << Dbl(pairForceT*100.0/time, 12 , 6, true) << std::endl;
      out << "Bond Forces          " 
          << Dbl(bondForceT*iStepInv, 12, 6) << " sec     " 
          << Dbl(bondForceT*ratio, 12, 6)    << " sec   " 
          << Dbl(bondForceT*100.0/time, 12 , 6, true) << std::endl;
      out << "Integrate2           " 
          << Dbl(integrate2T*iStepInv, 12, 6) << " sec     " 
          << Dbl(integrate2T*ratio, 12, 6)    << " sec   " 
          << Dbl(integrate2T*100.0/time, 12, 6, true) << std::endl;
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

   /*
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
