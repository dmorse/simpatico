/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PairEnergySelector.h"
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/accumulators/Average.h>                    // member template 
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   PairEnergySelector::PairEnergySelector(Simulation& simulation)
    : Analyzer(simulation),
      outputIntra_(),
      outputInter_(),
      pairPotentialPtr_(&simulation.pairPotential()),
      isInitialized_(false)
   {  setClassName("PairEnergySelector"); }

   /*
   * Destructor.
   */
   PairEnergySelector::~PairEnergySelector()
   { }

   /*
   * Load internal state from an archive.
   */
   void PairEnergySelector::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      nAtomType_ = simulation().nAtomType();

      selectors_.allocate(nAtomType_, 2 * nAtomType_);

      int type; int i; int j;
      int shift = 0;
      // Type = 1 -> INTRA ; Type = 2 -> INTER

      for (type = 1; type < 3; type++) {
         if (type == 2)
            shift = nAtomType_;

         // create a matrix of selectors INTRA and INTER side by side
         for (i = 0; i < nAtomType_; i++) {
            for (j = 0; j < nAtomType_; j++) {
               selectors_(i,j + shift).setAtom1TypeId(i);
               selectors_(i,j + shift).setAtom2TypeId(j);
               selectors_(i,j + shift).setPairType((PairSelector::PairType)type);
               selectors_(i,j + shift).setAvoidDoubleCounting(false);
            }
         }
      }

      ar >> pairEnergies_;

      if(simulation().domain().isMaster()) {
         // Two upper triangluar matrices
         simulation().fileMaster().
            openOutputFile(outputFileName("_intra.dat"), outputIntra_);
         simulation().fileMaster().
            openOutputFile(outputFileName("_inter.dat"), outputInter_);
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void PairEnergySelector::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << pairEnergies_;
   }

   /*
   * Read interval and outputFileName. 
   */
   void PairEnergySelector::readParameters(std::istream& in)
   {
//      int i = 0;
//      char hostname[256];
//      int myrank = simulation().domain().communicator().Get_rank();
//      gethostname(hostname, sizeof(hostname));
//      printf("rank %i PID %d on %s ready for attach\n", myrank, getpid(), hostname);
//      fflush(stdout);
//      sleep(15);

      readInterval(in);
      readOutputFileName(in);
      nAtomType_ = simulation().nAtomType();

      selectors_.allocate(nAtomType_, 2 * nAtomType_);

      int type; int i; int j;
      int shift = 0;
      // Type = 1 -> INTRA ; Type = 2 -> INTER

      for (type = 1; type < 3; type++) {
         if (type == 2)
            shift = nAtomType_;

         // create a matrix of selectors INTRA and INTER side by side
         for (i = 0; i < nAtomType_; i++) {
            for (j = 0; j < nAtomType_; j++) {
               selectors_(i,j + shift).setAtom1TypeId(i);
               selectors_(i,j + shift).setAtom2TypeId(j);
               selectors_(i,j + shift).setPairType((PairSelector::PairType)type);
               selectors_(i,j + shift).setAvoidDoubleCounting(false);
            }
         }
      }

      if(simulation().domain().isMaster()) {
         // Two upper triangluar matrices
         int nEnergies = ((nAtomType_ * (nAtomType_ - 1)) / 2) + nAtomType_;
         nEnergies *= 2;

         // add 2 (one for each sum)
         pairEnergies_.allocate(nEnergies + 2);

         simulation().fileMaster().
            openOutputFile(outputFileName("_intra.dat"), outputIntra_);
         simulation().fileMaster().
            openOutputFile(outputFileName("_inter.dat"), outputInter_);
      }
      isInitialized_ = true;
   }

   /*
   * Dump configuration to file
   */
   void PairEnergySelector::sample(long iStep)
   {
      if (isAtInterval(iStep))  {

         if (!pairPotentialPtr_->cellList().isValid()) {
            pairPotentialPtr_->buildCellList();
            pairPotentialPtr_->buildPairList();
         }
         pairListPtr_ = &pairPotentialPtr_->pairList();

         Vector f;
         double rsq;
         PairIterator iter;
         Atom*  atom1Ptr;
         Atom*  atom2Ptr;
         int    type0, type1;

         // Create energies array and initialize to zero
         DMatrix<double> energies;
         DMatrix<double> totalEnergies;
         energies.allocate(nAtomType_, 2 * nAtomType_);
         totalEnergies.allocate(nAtomType_, 2 * nAtomType_);

         // For each pair selector
         int i = 0; int j = 0;
         for (i = 0; i < nAtomType_; i++) {
            for (j = 0; j < nAtomType_; j++) {
               energies(i,j) = 0.0;
               totalEnergies(i,j) = 0.0;
               energies(i, j + nAtomType_) = 0.0;
               totalEnergies(i, j + nAtomType_) = 0.0;
            }
         }

         int myrank = simulation().domain().communicator().Get_rank();

         // For each pair in pairlist
         for (pairListPtr_->begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom1Ptr, atom2Ptr);

            for (i = 0; i < nAtomType_; i++) {
               for (j = 0; j < 2 * nAtomType_; j++) {
                  if (selectors_(i, j).match(*atom1Ptr, *atom2Ptr)) {
                     double rsq; double energy;

                     f.subtract(atom1Ptr->position(), atom2Ptr->position());
                     rsq = f.square();
                     energy = simulation().pairPotential().pairEnergy(rsq,
                           atom1Ptr->typeId(), atom2Ptr->typeId());

                     if (!atom2Ptr->isGhost()) {
                        energies(i,j) += energy;
                     }
                     else {
                        energies(i,j) += 0.5*energy;
                     }

                     goto nextatom;
                  }
               }
            }

            nextatom:rsq += 0;

         } // end pair loop

         // Collect energy matrices on master rank
         for (i = 0; i < nAtomType_; i++) {
            for (j = 0; j < 2 * nAtomType_; j++) {
               simulation().domain().communicator(). Reduce(&energies(i,j),
                     &totalEnergies(i,j), 1, MPI::DOUBLE, MPI::SUM, 0);
            }
         }

         if (simulation().domain().isMaster()) {
            int nEnergies = ((nAtomType_ * (nAtomType_ - 1)) / 2) + nAtomType_;
            nEnergies *= 2;
            int pairId = 0;
            double intraSum = 0.0;
            double interSum = 0.0;

            int idex = 0;
            for (i = 0; i < nAtomType_; i++) {
               for (j = 0; j < nAtomType_; j++) {
                  if (j >= i) {
                     double intraEnergy = totalEnergies(i, j);
                     double interEnergy = totalEnergies(i, j + nAtomType_);

                     // add the lower and upper parts of the matrices together
                     // so that INTER 1 2 and INTER 2 1 are reported as one.
                     if (i != j) {
                        intraEnergy += totalEnergies(j, i);
                        interEnergy += totalEnergies(j, i + nAtomType_);
                     }
                     pairEnergies_[idex].sample(intraEnergy);
                     pairEnergies_[idex + nEnergies/2].sample(interEnergy);
                     outputIntra_ << Dbl(intraEnergy);
                     outputInter_ << Dbl(interEnergy);
                     interSum += interEnergy;
                     intraSum += intraEnergy;

                     idex++;
                  }
               }
            }

            /* sample and output the sum */
            pairEnergies_[nEnergies].sample(intraSum);
            pairEnergies_[nEnergies+1].sample(interSum);
            outputIntra_ << Dbl(intraSum) << std::endl;
            outputInter_ << Dbl(interSum) << std::endl;
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void PairEnergySelector::output()
   {
      if (simulation().domain().isMaster()) {
         outputIntra_.close();
         outputInter_.close();

         int nEnergies = ((nAtomType_ * (nAtomType_ - 1)) / 2) + nAtomType_;
/*
         simulation().fileMaster().
               openOutputFile(outputFileName("_intra.prm"), outputIntra_);
         simulation().fileMaster().
               openOutputFile(outputFileName("_inter.prm"), outputInter_);
         ParamComposite::writeParam(outputIntra_);
         ParamComposite::writeParam(outputInter_);

         int pairId;
         for (pairId = 0; pairId < nEnergies; pairId++) {
            outputIntra_ << selectors_[pairId] << std::endl;
            outputInter_ << selectors_[pairId+nEnergies] << std::endl;
         }
         outputIntra_ << selectors_[2*nEnergies-1] << std::endl;
         outputInter_ << selectors_[2*nEnergies] << std::endl;

         outputIntra_.close();
         outputInter_.close();
*/
         simulation().fileMaster().
               openOutputFile(outputFileName("_intra.ave"), outputIntra_);
         simulation().fileMaster().
               openOutputFile(outputFileName("_inter.ave"), outputInter_);

         int idex = 0;
         int i; int j;
         for (i = 0; i < nAtomType_; i++) {
            for (j = 0; j < nAtomType_; j++) {
               if (j >= i) {
                  outputIntra_ << "INTRA  " << i << "  " << j << std::endl;
                  outputInter_ << "INTER  " << i << "  " << j << std::endl;
                  pairEnergies_[idex].output(outputIntra_);
                  pairEnergies_[idex + nEnergies].output(outputInter_);

                  idex++;
               }
            }
         }
         outputIntra_ << "SUM" << std::endl;
         outputInter_ << "SUM" << std::endl;
         pairEnergies_[2*nEnergies].output(outputIntra_);
         pairEnergies_[2*nEnergies+1].output(outputInter_);

         outputIntra_.close();
         outputInter_.close();
      }
   }
}
