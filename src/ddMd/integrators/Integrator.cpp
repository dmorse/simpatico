/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Integrator.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GroupStorage.tpp>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/analyzers/AnalyzerManager.h>
#include <ddMd/potentials/pair/PairPotential.h>
#ifdef SIMP_BOND
#include <ddMd/potentials/bond/BondPotential.h>
#endif
#ifdef SIMP_ANGLE
#include <ddMd/potentials/angle/AnglePotential.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef SIMP_EXTERNAL
#include <ddMd/potentials/external/ExternalPotential.h>
#endif

#include <util/ensembles/BoundaryEnsemble.h>
#include <util/mpi/MpiLoader.h>
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
       isSetup_(false),
       saveFileName_(),
       saveInterval_(0)
   {}

   /*
   * Destructor.
   */
   Integrator::~Integrator()
   {}

   /*
   * Read saveInterval and saveFileName.
   */
   void Integrator::readParameters(std::istream& in)
   {
      read<int>(in, "saveInterval", saveInterval_);
      if (saveInterval_ > 0) {
         if (Analyzer::baseInterval > 0) {
            if (saveInterval_ % Analyzer::baseInterval != 0) {
               UTIL_THROW("saveInterval is not a multiple of baseInterval");
            }
         } else {
            UTIL_THROW("Analyzer::baseInterval is not positive");
         }
         read<std::string>(in, "saveFileName", saveFileName_);
      }
   }

   /*
   * Load saveInterval and saveFileName from restart archive.
   */
   void Integrator::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<int>(ar, "saveInterval", saveInterval_);
      if (saveInterval_ > 0) {
         if (Analyzer::baseInterval > 0) {
            if (saveInterval_ % Analyzer::baseInterval != 0) {
               UTIL_THROW("saveInterval is not a multiple of baseInterval");
            }
         } else {
            UTIL_THROW("Analyzer::baseInterval is not positive");
         }
         loadParameter<std::string>(ar, "saveFileName", saveFileName_);
      }

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(iStep_);
      loader.load(isSetup_);
   }

   /*
   * Save saveInterval and saveFileName to restart archive.
   */
   void Integrator::save(Serializable::OArchive& ar)
   {
      ar << saveInterval_;
      if (saveInterval_ > 0) {
         ar << saveFileName_;
      }
      ar << iStep_;
      ar << isSetup_;
   }

   /*
   * Exchange atoms, build PairList and compute forces.
   */
   void Integrator::setupAtoms()
   {
      // Precondition
      if (atomStorage().isCartesian()) {
         UTIL_THROW("Atom coordinates are Cartesian");
      }

      atomStorage().clearSnapshot();
      exchanger().exchange();
      pairPotential().buildCellList();
      atomStorage().transformGenToCart(boundary());
      atomStorage().makeSnapshot();
      pairPotential().buildPairList();
      if (simulation().boundaryEnsemble().isRigid()) {
         simulation().computeForces();
      } else {
         simulation().computeForcesAndVirial();
      }

      // Postcondition - coordinates are Cartesian
      if (!atomStorage().isCartesian()) {
         UTIL_THROW("Atom coordinates are not Cartesian");
      }
   }

   /*
   * Compute forces for all atoms, with timing.
   */
   void Integrator::computeForces()
   {
      // Precondition
      if (!atomStorage().isCartesian()) {
         UTIL_THROW("Atom coordinates are not Cartesian");
      }

      timer_.stamp(MISC);
      simulation().zeroForces();
      timer_.stamp(ZERO_FORCE);
      pairPotential().computeForces();
      timer_.stamp(PAIR_FORCE);
      #ifdef SIMP_BOND
      if (nBondType()) {
         bondPotential().computeForces();
         timer_.stamp(BOND_FORCE);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType()) {
         anglePotential().computeForces();
         timer_.stamp(ANGLE_FORCE);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType()) {
         dihedralPotential().computeForces();
         timer_.stamp(DIHEDRAL_FORCE);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal()) {
         externalPotential().computeForces();
         timer_.stamp(EXTERNAL_FORCE);
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
      // Precondition
      if (!atomStorage().isCartesian()) {
         UTIL_THROW("Atom coordinates are not Cartesian");
      }

      timer_.stamp(MISC);
      simulation().zeroForces();
      timer_.stamp(ZERO_FORCE);
      pairPotential().computeForcesAndStress(domain().communicator());
      timer_.stamp(PAIR_FORCE);
      #ifdef SIMP_BOND
      if (nBondType()) {
         bondPotential().computeForcesAndStress(domain().communicator());
         timer_.stamp(BOND_FORCE);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType()) {
         anglePotential().computeForcesAndStress(domain().communicator());
         timer_.stamp(ANGLE_FORCE);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType()) {
         dihedralPotential().computeForcesAndStress(domain().communicator());
         timer_.stamp(DIHEDRAL_FORCE);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal()) {
         externalPotential().computeForces();
         timer_.stamp(EXTERNAL_FORCE);
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
      double factor2 = double(nProc)/(double(iStep_)*double(nAtomTot));
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
      double transformFT = timer().time(TRANSFORM_F);
      totalT += transformFT;
      out << "Transform (forward)  " 
          << Dbl(transformFT*factor1, 12, 6)
          << "   "
          << Dbl(transformFT*factor2, 12, 6)
          << "   " << Dbl(100.0*transformFT/time, 12, 6, true) << std::endl;
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
      double transformRT = timer().time(TRANSFORM_R);
      totalT += transformRT;
      out << "Transform (reverse)  " 
          << Dbl(transformRT*factor1, 12, 6)
          << "   "
          << Dbl(transformRT*factor2, 12, 6)
          << "   " << Dbl(100.0*transformRT/time, 12, 6, true) << std::endl;
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
      double zeroForceT = timer().time(ZERO_FORCE);
      totalT += zeroForceT;
      out << "Zero Forces          " 
          << Dbl(zeroForceT*factor1, 12, 6)
          << "   "
          << Dbl(zeroForceT*factor2, 12, 6)
          << "   " << Dbl(100.0*zeroForceT/time, 12 , 6, true) << std::endl;
      double pairForceT = timer().time(PAIR_FORCE);
      totalT += pairForceT;
      out << "Pair Forces          " 
          << Dbl(pairForceT*factor1, 12, 6)
          << "   "
          << Dbl(pairForceT*factor2, 12, 6)
          << "   " << Dbl(100.0*pairForceT/time, 12 , 6, true) << std::endl;
      #ifdef SIMP_BOND
      if (nBondType()) {
         double bondForceT = timer().time(BOND_FORCE);
         totalT += bondForceT;
         out << "Bond Forces          " 
             << Dbl(bondForceT*factor1, 12, 6)
             << "   "
             << Dbl(bondForceT*factor2, 12, 6)
             << "   " << Dbl(100.0*bondForceT/time, 12 , 6, true) << std::endl;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType()) {
         double angleForceT = timer().time(ANGLE_FORCE);
         totalT += angleForceT;
         out << "Angle Forces         " 
             << Dbl(angleForceT*factor1, 12, 6)
             << "   "
             << Dbl(angleForceT*factor2, 12, 6)
             << "   " << Dbl(100.0*angleForceT/time, 12 , 6, true) << std::endl;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType()) {
         double dihedralForceT = timer().time(DIHEDRAL_FORCE);
         totalT += dihedralForceT;
         out << "Dihedral Forces      " 
             << Dbl(dihedralForceT*factor1, 12, 6) 
             << "   "
             << Dbl(dihedralForceT*factor2, 12, 6)
             << "   " << Dbl(100.0*dihedralForceT/time, 12 , 6, true) << std::endl;
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal()) {
         double externalForceT = timer().time(DIHEDRAL_FORCE);
         totalT += externalForceT;
         out << "External Forces      " 
             << Dbl(externalForceT*factor1, 12, 6) 
             << "   "
             << Dbl(externalForceT*factor2, 12, 6)
             << "   " << Dbl(100.0*externalForceT/time, 12 , 6, true) << std::endl;
      }
      #endif
      double integrate2T = timer().time(INTEGRATE2);
      totalT += integrate2T;
      out << "Integrate2           " 
          << Dbl(integrate2T*factor1, 12, 6) 
          << "   "
          << Dbl(integrate2T*factor2, 12, 6) 
          << "   " << Dbl(100.0*integrate2T/time, 12, 6, true) << std::endl;
      double analyzerT = timer().time(ANALYZER);
      totalT += analyzerT;
      out << "Analyzers            " 
          << Dbl(analyzerT*factor1, 12, 6)
          << "   "
          << Dbl(analyzerT*factor2, 12, 6)
          << "   " << Dbl(100.0*analyzerT/time, 12, 6, true) << std::endl;
      #ifdef DDMD_MODIFIERS
      double modifierT = timer().time(MODIFIER);
      totalT += modifierT;
      out << "Modifiers            " 
          << Dbl(modifierT*factor1, 12, 6)
          << "   "
          << Dbl(modifierT*factor2, 12, 6)
          << "   " << Dbl(100.0*modifierT/time, 12, 6, true) << std::endl;
      #endif
      out << std::endl;

      // Output info about timer resolution
      double tick = MPI::Wtick();
      out << "Timer resolution     " 
          << Dbl(tick, 12, 6) 
          << "   "
          << Dbl(tick*double(nProc)/double(nAtomTot), 12, 6) 
          << "   " << Dbl(100.0*tick*double(iStep_)/time, 12, 6, true) 
          << std::endl;
      out << std::endl;

      // Output exchange / reneighbor statistics
      int buildCounter = pairPotential().pairList().buildCounter(); 
      out << "buildCounter         " 
                  << Int(buildCounter, 12)
                  << std::endl;
      out << "steps per build      "
                  << Dbl(double(iStep_)/double(buildCounter), 12, 6)
                  << std::endl;
      out << std::endl;

   }

   /*
   * Clear timing, dynamical state, statistics, and analyzer accumulators.
   */
   void Integrator::clear()
   { 
      iStep_ = 0;
      initDynamicalState();
      timer().clear(); 
      simulation().exchanger().timer().clear();
      simulation().buffer().clearStatistics();
      atomStorage().clearStatistics();
      pairPotential().pairList().clearStatistics();
      #ifdef SIMP_BOND
      bondStorage().clearStatistics();
      #endif
      #ifdef SIMP_ANGLE
      angleStorage().clearStatistics();
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralStorage().clearStatistics();
      #endif
      simulation().analyzerManager().clear();
   }

}
