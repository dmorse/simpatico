/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbReptateMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcSimulation/mc_potentials.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/chemistry/Bond.h>

#include <simp/species/Linear.h>

#include <util/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor
   */
   CfbReptateMove::CfbReptateMove(McSystem& system) :
      CfbLinear(system),
      bondTypeId_(-1),
      nJunction_(0),
      junctions_(),
      lTypes_(),
      uTypes_(),
      maskPolicy_(MaskBonded),
      hasAutoCorr_(false),
      autoCorrCapacity_(0),
      outputFileName_(),
      accumulators_()
   {  setClassName("CfbReptateMove"); }

   /*
   * Read parameters and initialize.
   */
   void CfbReptateMove::readParameters(std::istream& in)
   {
      // Read parameters
      readProbability(in);
      CfbLinear::readParameters(in); // read speciesId_ and nTrial_

      hasAutoCorr_ = false; // default value
      readOptional<bool>(in, "hasAutoCorr", hasAutoCorr_);
      if (hasAutoCorr_) {
         read<int>(in, "autoCorrCapacity", autoCorrCapacity_);
         read<std::string>(in, "outputFileName", outputFileName_);
      }

      // Identify bond type, check that it is same for all bonds
      Species* speciesPtr = &simulation().species(speciesId());
      int nBond = speciesPtr->nBond();
      bondTypeId_ = speciesPtr->speciesBond(0).typeId();
      for (int i = 1; i < nBond; ++i) {
         if (bondTypeId_ != speciesPtr->speciesBond(i).typeId()) {
            UTIL_THROW("Unequal bond type ids");
         }
      }

      #ifndef SIMP_NOPAIR
      // Identify policy for masking nonbonded interactions.
      maskPolicy_ = simulation().maskedPairPolicy();
      #endif

      // Allocate memory for junction arrays
      int nAtom = speciesPtr->nAtom();
      junctions_.allocate(nAtom-1);
      lTypes_.allocate(nAtom-1);
      uTypes_.allocate(nAtom-1);

      // Identify junctions between atoms of different types
      int lType, uType;
      nJunction_ = 0;
      for (int i = 0; i < nAtom - 1; ++i) {
         lType = speciesPtr->atomTypeId(i);
         uType = speciesPtr->atomTypeId(i+1);
         if (lType != uType) {
            junctions_[nJunction_] = i;
            lTypes_[nJunction_] = lType;
            uTypes_[nJunction_] = uType;
            ++nJunction_;
         }
      }

      if (hasAutoCorr_) {
         // Allocate memory for autocorrelation accumulators
         int moleculeCapacity = simulation().species(speciesId()).capacity();
         accumulators_.allocate(moleculeCapacity);
         for (int i = 0; i < moleculeCapacity; i++) {
            accumulators_[i].setParam(autoCorrCapacity_);
         }
      }

   }

   /*
   * Load from archive.
   */
   void CfbReptateMove::loadParameters(Serializable::IArchive& ar)
   {
      // Read parameters
      McMove::loadParameters(ar);
      CfbLinear::loadParameters(ar);

      // Read autocorrelation parameters
      loadParameter<bool>(ar, "hasAutoCorr", hasAutoCorr_);
      if (hasAutoCorr_) {
         loadParameter<int>(ar, "autoCorrCapacity", autoCorrCapacity_);
         loadParameter<std::string>(ar, "outputFileName", outputFileName_);
         ar & accumulators_;
      }
      ar & bondTypeId_;
      ar & maskPolicy_;
      ar & nJunction_;
      ar & junctions_;
      ar & lTypes_;
      ar & uTypes_;

      // Validate
      Species* speciesPtr = &simulation().species(speciesId());
      int nBond = speciesPtr->nBond();
      for (int i = 0; i < nBond; ++i) {
         if (bondTypeId_ != speciesPtr->speciesBond(i).typeId()) {
            UTIL_THROW("Inconsistent or unequal bond type ids");
         }
      }
      #ifndef SIMP_NOPAIR
      if (maskPolicy_ != simulation().maskedPairPolicy()) {
         UTIL_THROW("Inconsistent values of maskPolicy_");
      }
      #endif
   }

   /*
   * Save to archive.
   */
   void CfbReptateMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      CfbLinear::save(ar);
      ar & hasAutoCorr_;
      if (hasAutoCorr_) {
         ar & autoCorrCapacity_;
         ar & outputFileName_;
         ar & accumulators_;
      }
      ar & bondTypeId_;
      ar & maskPolicy_;
      ar & nJunction_;
      ar & junctions_;
      ar & lTypes_;
      ar & uTypes_;
   }

   /*
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool CfbReptateMove::move()
   {
      Vector oldPos, newPos;
      double rosen_r,  rosen_f, ratio;
      double energy_r, energy_f;
      Atom *tailPtr; // pointer to the tail atom (to be removed)
      Atom *atomPtr; // resetable atom pointer
      Molecule *molPtr; // pointer to chosen molecule
      int nAtom, sign, headId, tailId, oldType, i;
      bool accept;

      incrementNAttempt();

      // Choose a molecule of specific species at random
      molPtr = &(system().randomMolecule(speciesId()));
      nAtom = molPtr->nAtom();

      // Choose which chain end to regrow.
      // "tail" = deletion end, "head" = addition end
      if (random().uniform(0.0, 1.0) > 0.5) {
         sign = +1;
         headId = nAtom - 1;
         tailId = 0;
      } else {
         sign = -1;
         headId = 0;
         tailId = nAtom - 1;
      }

      // Store current position and type of tail atom
      tailPtr = &(molPtr->atom(tailId));
      oldPos = tailPtr->position();
      oldType = tailPtr->typeId();

      // Delete tail atom
      deleteAtom(*molPtr, tailId, -sign, rosen_r, energy_r);
      #ifndef SIMP_NOPAIR
      // Delete from McSystem cell list
      system().pairPotential().deleteAtom(*tailPtr);
      #endif

      // Regrow head, using tailPtr to store properties of the new head.
      atomPtr = &(molPtr->atom(headId)); // new pivot atom
      tailPtr->setTypeId(atomPtr->typeId()); // new head atom
      tailPtr->mask().clear();
      if (maskPolicy_ == MaskBonded) {
         tailPtr->mask().append(*atomPtr);
      }
      addAtom(*molPtr, *tailPtr, *atomPtr, headId, sign, rosen_f, energy_f);

      // Calculate junction factor for heteropolymers.
      double jFactor;
      if (nJunction_ > 0) {
         jFactor = junctionFactor(molPtr, sign);
         // Upon return, type ids are restored to original values.
      } else {
         jFactor = 1.0;
      }

      // Decide whether to accept or reject
      ratio = jFactor * rosen_f / rosen_r;
      accept = random().metropolis(ratio);

      // Restore original type and connectivity mask for tail Atom
      tailPtr->setTypeId(oldType);
      tailPtr->mask().clear();
      if (maskPolicy_ == MaskBonded) {
         atomPtr = tailPtr + sign;
         tailPtr->mask().append(*atomPtr);
      }

      if (accept) {

         // Store position of new head atom
         newPos = tailPtr->position();

         // Shift position of tail
         atomPtr = tailPtr + sign;
         tailPtr->position() = atomPtr->position();
         #ifndef SIMP_NOPAIR
         // Add tail back to system cell list
         system().pairPotential().addAtom(*tailPtr);
         #endif

         // Shift atom positions towards the head
         for (i = 1; i < nAtom - 1; ++i) {
            atomPtr->position() = (atomPtr+sign)->position();
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*atomPtr);
            #endif
            atomPtr += sign;
         }
         assert(atomPtr == &molPtr->atom(headId));

         // Move head atom to new chosen position
         atomPtr->position() = newPos;
         #ifndef SIMP_NOPAIR
         system().pairPotential().updateAtomCell(*atomPtr);
         #endif

         // Increment the number of accepted moves.
         incrementNAccept();

         if (hasAutoCorr_) {
            int molId = molPtr->id();
            double rsign = (double) sign;
            accumulators_[molId].sample(rsign);
            //accumulators_[molPtr->id()].sample((double) sign);
         }

      } else {

         // Restore old position of tail
         tailPtr->position() = oldPos;
         #ifndef SIMP_NOPAIR
         // Add tail back to System cell list.
         system().pairPotential().addAtom(*tailPtr);
         #endif

      }

      return accept;
   }

   /**
   * Calculate Boltzmann factor associated with all junctions.
   */
   double CfbReptateMove::junctionFactor(Molecule* molPtr, int sign)
   {
      double oldEnergy, newEnergy, factor;
      Atom* hAtomPtr; // Pointer to atom nearer head
      int hType;      // type id of atom nearer head
      int tType;      // type id of atom nearer tail
      int i;          // junction index
      int j;          // id of atom at below junction (lower)

      // Calculate factor by looping over junctions
      factor = 1.0;
      for (i = 0; i < nJunction_; ++i) {
         j = junctions_[i];
         if (sign == 1) {
            hAtomPtr = &(molPtr->atom(j+1));
            tType = lTypes_[i];
         } else {
            hAtomPtr = &(molPtr->atom(j));
            tType = uTypes_[i];
         }

         #ifndef SIMP_NOPAIR
         oldEnergy = system().pairPotential().atomEnergy(*hAtomPtr);
         #else
         oldEnergy = 0.0;
         #endif
         #ifdef SIMP_EXTERNAL
         if (system().hasExternalPotential()) {
            oldEnergy += system().externalPotential().atomEnergy(*hAtomPtr);
         }
         #endif

         #ifndef SIMP_NOPAIR
         hAtomPtr->setTypeId(tType);
         newEnergy = system().pairPotential().atomEnergy(*hAtomPtr);
         #else
         newEnergy = 0.0;
         #endif
         #ifdef SIMP_EXTERNAL
         if (system().hasExternalPotential()) {
            newEnergy += system().externalPotential().atomEnergy(*hAtomPtr);
         }
         #endif

         factor *= boltzmann(newEnergy - oldEnergy);
      }

      // Revert modified atom type Ids to original values
      #ifndef SIMP_NOPAIR
      for (i = 0; i < nJunction_; ++i) {
         j = junctions_[i];
         if (sign == 1) {
            hAtomPtr = &(molPtr->atom(j+1));
            hType    = uTypes_[i];
         } else {
            hAtomPtr = &(molPtr->atom(j));
            hType    = lTypes_[i];
         }
         hAtomPtr->setTypeId(hType);
      }
      return factor;

      #else //ifdef SIMP_NOPAIR
      return 1.0;
      #endif
   }

   /**
   * Output statistics about accepted reptation steps
   */
   void CfbReptateMove::output()
   {
      if (hasAutoCorr_) {
         DArray< double > autoCorrAvg;
         autoCorrAvg.allocate(autoCorrCapacity_);

         // Compute average autocorrelation function over molecules
         for (int i = 0; i < autoCorrCapacity_; i++) {
            int nAvg = 0;
            autoCorrAvg[i] = 0;
            for (int j = 0; j < system().nMolecule(speciesId()); j++) {
               if (accumulators_[j].nSample() > i) {
                  nAvg++;
                  autoCorrAvg[i] += accumulators_[j].autoCorrelation(i);
               }
            }
            autoCorrAvg[i] /= nAvg;
         }

         // Write out average autocorrelation
         std::ofstream outputFile;
         std::string fileName = outputFileName_;
         fileName += ".dat";
         system().fileMaster().openOutputFile(fileName, outputFile);
         for (int i = 0; i <  autoCorrCapacity_; i++) {
            outputFile << Int(i);
            write<double>(outputFile, autoCorrAvg[i]);
            outputFile << std::endl;
         }
         outputFile.close();

         // Sum over autocorrelation (only positive time lags)
         // Count first element with a factor 1/2, due to symmetry
         double acSum = 0;
         acSum = 0.5*autoCorrAvg[0];
         for (int i = 1; i < autoCorrCapacity_; i++) {
            acSum += autoCorrAvg[i];
         }

         // Write sum in .ave file
         fileName = outputFileName_;
         fileName += ".ave";
         system().fileMaster().openOutputFile(fileName, outputFile);
         outputFile << "autoCorrSum " << acSum << std::endl;
         outputFile.close();
      }
   }

}
