#ifndef  INTER_NOPAIR
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdPairEnergySelector.h"

#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/neighbor/PairList.h>
#include <mcMd/neighbor/PairIterator.h>
#include <mcMd/species/Species.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/archives/Serializable_includes.h>
#include <util/misc/FileMaster.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdPairEnergySelector::MdPairEnergySelector(MdSystem& system)
    : SystemAnalyzer<MdSystem>(system),
    pairListPtr_(&system.pairPotential().pairList()),
    pairPotentialPtr_(&system.pairPotential()),
    boundaryPtr_(&system.boundary())
   {  setClassName("MdPairEnergySelector"); }

   /*
   * Destructor
   */
   MdPairEnergySelector::~MdPairEnergySelector()
   {}


   /*
   * Read input parameters.
   */
   void MdPairEnergySelector::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "nPairTypes", nPairTypes_);
      selectors_.allocate(nPairTypes_);

      /* allocate one extra for sum */
      pairEnergies_.allocate(nPairTypes_ + 1);

      readDArray<PairSelector>(in, "selectors", selectors_, nPairTypes_);

      int pairId;
      for (pairId = 0; pairId < nPairTypes_; pairId++) {
         selectors_[pairId].setAvoidDoubleCounting(false);
      }

      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void MdPairEnergySelector::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);

      ar & nPairTypes_;
      selectors_.allocate(nPairTypes_);
      pairEnergies_.allocate(nPairTypes_ + 1);

      int i;
      for (i = 0; i < nPairTypes_; i++) {
         ar & selectors_[i];
      }
      for (i = 0; i <= nPairTypes_; i++) {
         ar & pairEnergies_[i];
      }


//      ar & selectors_;
//      ar & pairEnergies_;

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void MdPairEnergySelector::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Evaluate contributions to accumulators.
   */
   void MdPairEnergySelector::sample(long iStep)
   {
      if (isAtInterval(iStep)) {
         #ifndef INTER_NOPAIR

            // Loop over pairs and construct neighbor list per molecule
            PairIterator iter;
            Atom *atom1Ptr;
            Atom *atom2Ptr;

            // Create energies array and initialize to zero
            double energies[nPairTypes_];
            int pairId;
            for (pairId = 0; pairId < nPairTypes_; pairId++) {
               energies[pairId] = 0.0;
            }

            for (pairListPtr_->begin(iter); iter.notEnd(); ++iter) {
               iter.getPair(atom1Ptr, atom2Ptr);

               for (pairId = 0; pairId < nPairTypes_; pairId++) {
                  if (selectors_[pairId].match(*atom1Ptr, *atom2Ptr)) {
                     double rsq; double energy;

                     rsq = boundaryPtr_->distanceSq(atom1Ptr->position(),
                           atom2Ptr->position());
                     energy = pairPotentialPtr_->energy(rsq,
                           atom1Ptr->typeId(), atom2Ptr->typeId());

                     energies[pairId] += energy;
                  }
               }

            } // end pair loop

            double pairSum = 0.0;
            for (pairId = 0; pairId < nPairTypes_; pairId++) {
               pairEnergies_[pairId].sample(energies[pairId]);
               outputFile_ << Dbl(energies[pairId]);
               pairSum += energies[pairId];
            }

            /* sample the sum */
            pairEnergies_[nPairTypes_].sample(pairSum);
            outputFile_ << Dbl(pairSum);

            outputFile_ << std::endl;

         #endif
//         outputFile_ << std::endl;
      }
   }
 
   /* 
   * Summary
   */
   void MdPairEnergySelector::output() 
   {
      int pairId;
      // Close *.dat file
      outputFile_.close();

      // Open and write summary file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_ << std::endl;

      outputFile_ << "File format:" << std::endl;
      #ifndef INTER_NOPAIR
//      outputFile_ << "  ";
      for (pairId = 0; pairId < nPairTypes_; pairId++) {
         outputFile_ << selectors_[pairId] << std::endl;
      }
      #endif
      outputFile_ << "SUM" << std::endl;
      outputFile_ << std::endl;

      for (pairId = 0; pairId <= nPairTypes_; pairId++) {
         pairEnergies_[pairId].output(outputFile_);
         outputFile_ << std::endl;
      }


      outputFile_.close();

   }

}
#endif 
