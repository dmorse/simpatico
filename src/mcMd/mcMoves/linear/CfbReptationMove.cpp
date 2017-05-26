/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbReptationMove.h"
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

   /* 
   * Constructor
   */
   CfbReptationMove::CfbReptationMove(McSystem& system) : 
      CfbEndBase(system),
      speciesId_(-1),
      bondTypeId_(-1),
      nJunction_(0),
      junctions_(),
      lTypes_(),
      uTypes_(),
      maskPolicy_(MaskBonded),
      hasAutoCorr_(0),
      autoCorrCapacity_(0),
      outputFileName_(),
      acceptedStepsAccumulators_()
   {  setClassName("CfbReptationMove"); } 
   
   /* 
   * Read parameters speciesId, nRegrow, and nTrial
   * also read hasAutoCorr, which may be either 0 or 1
   * if 1, also read autoCorrCapacity and autoCorrName
   */
   void CfbReptationMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nTrial", nTrial_);

      // Check that 0 < nTrial_ <= MaxTrial
      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }

      // Check that the Species is a subclass of Linear 
      Linear* speciesPtr;
      speciesPtr = dynamic_cast<Linear*>(&(simulation().species(speciesId_)));
      if (!speciesPtr) {
         UTIL_THROW("Species is not a subclass of Linear");
      }

      // Identify bond type 
      bondTypeId_ = speciesPtr->speciesBond(0).typeId();

      // Check that bond type is the same for all bonds
      int nBond = speciesPtr->nBond();
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
            //std::cout << nJunction_ << "  " << i << "  " 
                     // << lType << "  " << uType << std::endl;
            ++nJunction_;
         }
      }

      read<int>(in, "hasAutoCorr", hasAutoCorr_);
 
      // Read autocorrelation parameters 
      if (hasAutoCorr_) {
         read<int>(in, "autoCorrCapacity", autoCorrCapacity_);
         read<std::string>(in, "autoCorrName", outputFileName_);

         // Allocate memory for autocorrelation function of reptation steps
         int moleculeCapacity = simulation().species(speciesId_).capacity();
         acceptedStepsAccumulators_.allocate(moleculeCapacity); 
         for (int i = 0; i < moleculeCapacity; i++)
         {
            acceptedStepsAccumulators_[i].setParam(autoCorrCapacity_);
         }
      }

   }

   /* 
   * Load from archive.
   */
   void CfbReptationMove::loadParameters(Serializable::IArchive& ar) 
   {
      // Read parameters
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "nTrial", nTrial_);
      ar & bondTypeId_;
      ar & maskPolicy_;

      // Validate
      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }
      Linear* speciesPtr;
      speciesPtr = dynamic_cast<Linear*>(&(simulation().species(speciesId_)));
      if (!speciesPtr) {
         UTIL_THROW("Species is not a subclass of Linear");
      }
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

      ar & junctions_; 
      ar & lTypes_; 
      ar & uTypes_; 

      // Read autocorrelation parameters 
      loadParameter<int>(ar, "hasAutoCorr", hasAutoCorr_);
      if (hasAutoCorr_) {
         loadParameter<int>(ar, "autoCorrCapacity", autoCorrCapacity_);
         loadParameter<std::string>(ar, "autoCorrName", outputFileName_);
         ar & acceptedStepsAccumulators_;
      }

   }

   /* 
   * Load from archive.
   */
   void CfbReptationMove::save(Serializable::OArchive& ar) 
   {
      McMove::save(ar);
      ar & speciesId_;
      ar & nTrial_;
      ar & bondTypeId_;
      ar & maskPolicy_;
      ar & junctions_; 
      ar & lTypes_; 
      ar & uTypes_; 
      ar & hasAutoCorr_;
      if (hasAutoCorr_) {
         ar & autoCorrCapacity_;
         ar & outputFileName_;
         ar & acceptedStepsAccumulators_;
      }
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool CfbReptationMove::move() 
   {
      Vector    oldPos, newPos;
      double    rosen_r,  rosen_f, ratio;
      double    energy_r, energy_f;
      Atom     *tailPtr; // pointer to the tail atom (to be removed)
      Atom     *atomPtr; // resetable atom pointer
      Molecule *molPtr;  // pointer to randomly chosen molecule
      int       length, sign, headId, tailId, oldType, i;
      bool      accept;
     
      incrementNAttempt();

      // Choose a molecule at random
      molPtr = &(system().randomMolecule(speciesId_));
      length = molPtr->nAtom();

      // Choose which chain end to regrow
      if (random().uniform(0.0, 1.0) > 0.5) {
         sign = +1;
         headId = length - 1;
         tailId = 0;
      } else {
         sign = -1;
         headId = 0;
         tailId = length - 1;
      }

      // Store current position and type of tail atom
      tailPtr = &(molPtr->atom(tailId));
      oldPos  = tailPtr->position();
      oldType = tailPtr->typeId();
   
      // Delete tail atom
      atomPtr = tailPtr + sign;
      deleteEndAtom(tailPtr, atomPtr, bondTypeId_, rosen_r, energy_r);
      #ifndef SIMP_NOPAIR
      // Delete from McSystem cell list
      system().pairPotential().deleteAtom(*tailPtr);
      #endif
   
      // Regrow head, using tailPtr to store properties of the new head.
      atomPtr = &(molPtr->atom(headId));
      tailPtr->setTypeId(atomPtr->typeId());
      tailPtr->mask().clear();
      if (maskPolicy_ == MaskBonded) {
         tailPtr->mask().append(*atomPtr);
      }
      addEndAtom(tailPtr, atomPtr, bondTypeId_, rosen_f, energy_f);

      // Calculate junction factor for heteropolymers.
      double jFactor = 1.0;
      if (nJunction_ > 0) {
         jFactor = junctionFactor(molPtr, sign);
         // Upon return, type ids are restored to original values.
      }

      // Decide whether to accept or reject
      ratio  = jFactor * rosen_f / rosen_r;
      accept = random().metropolis(ratio);

      // Restore original type and connectivity mask for tail Atom
      tailPtr->setTypeId(oldType);
      tailPtr->mask().clear();
      if (maskPolicy_ == MaskBonded) {
         atomPtr = tailPtr + sign;
         tailPtr->mask().append(*atomPtr);
      }

      if (accept) {

         // Store new head position
         newPos = tailPtr->position();
  
         // Shift position of tail 
         atomPtr = tailPtr + sign;
         tailPtr->position() = atomPtr->position();

         #ifndef SIMP_NOPAIR
         // Add tail back to system cell list
         system().pairPotential().addAtom(*tailPtr);
         #endif
   
         // Shift atom positions towards the head
         for (i=1; i < length - 1; ++i) {
            //system().moveAtom(*atomPtr, (atomPtr+sign)->position());
            atomPtr->position() = (atomPtr+sign)->position();
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*atomPtr);
            #endif
            atomPtr += sign;
         }
   
         // Move head atom to new chosen position
         // system().moveAtom(*atomPtr, newPos);
         atomPtr->position() = newPos;
         #ifndef SIMP_NOPAIR
         system().pairPotential().updateAtomCell(*atomPtr);
         #endif

         // Increment the number of accepted moves.
         incrementNAccept();

         if (hasAutoCorr_)
         {
            // Store step direction in autocorrelator
            acceptedStepsAccumulators_[molPtr->id()].sample((double) sign);
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
   double CfbReptationMove::junctionFactor(Molecule* molPtr, int sign) 
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
            tType    = lTypes_[i];
         } else {
            hAtomPtr = &(molPtr->atom(j));
            tType    = uTypes_[i];
         }

         #ifndef SIMP_NOPAIR
         oldEnergy = system().pairPotential().atomEnergy(*hAtomPtr);
         #endif

         #ifdef SIMP_EXTERNAL
         if (system().hasExternalPotential()) {
            oldEnergy += system().externalPotential().atomEnergy(*hAtomPtr);
         }
         #endif

         #ifndef SIMP_NOPAIR
         hAtomPtr->setTypeId(tType);
         newEnergy = system().pairPotential().atomEnergy(*hAtomPtr);
         #endif

         #ifdef SIMP_EXTERNAL
         if (system().hasExternalPotential()) {
            newEnergy += system().externalPotential().atomEnergy(*hAtomPtr);
         }
         #endif

         #ifndef SIMP_NOPAIR
         factor *= boltzmann(newEnergy - oldEnergy);
         #endif
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
   void CfbReptationMove::output()
   {
      if (hasAutoCorr_)
      {
         DArray< double > autoCorrAvg;
         autoCorrAvg.allocate(autoCorrCapacity_);

         // Calculate average autocorrelation function over all
         // accumulators
         for (int i = 0; i < autoCorrCapacity_; i++) {
            int nAvg = 0;
            autoCorrAvg[i] = 0;
            for (int j = 0; j < system().nMolecule(speciesId_); j++)
            {
               if (acceptedStepsAccumulators_[j].nSample() > i)
               {
                  nAvg++;
                  autoCorrAvg[i] += acceptedStepsAccumulators_[j].
                     autoCorrelation(i);
               }
            }
            autoCorrAvg[i] /= nAvg;
         }

         std::ofstream outputFile;

         // Write out average autocorrelation
         system().fileMaster().openOutputFile(outputFileName_+".dat",
            outputFile);

         for (int i = 0; i <  autoCorrCapacity_; i++)
         {
            outputFile << Int(i);
            write<double>(outputFile, autoCorrAvg[i]);
            outputFile << std::endl;
         }
      
         outputFile.close();

         // Sum over autocorrelation (only positive time lags)
         // Count first element with a factor 1/2, due to symmetry
         double acSum = 0;
         acSum = 0.5*autoCorrAvg[0];
         for (int i = 1; i < autoCorrCapacity_; i++)
            acSum += autoCorrAvg[i];

         // write sum in .ave file
         system().fileMaster().openOutputFile(outputFileName_+".ave",
            outputFile);
         outputFile << "autoCorrSum " << acSum << std::endl;
         outputFile.close();
      }
   }

}
