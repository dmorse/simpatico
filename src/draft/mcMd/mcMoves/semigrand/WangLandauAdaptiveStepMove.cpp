/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WangLandauAdaptiveStepMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/species/GeneralpolymerSG.h>

#include <mcMd/chemistry/Molecule.h> //testing
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   WangLandauAdaptiveStepMove::WangLandauAdaptiveStepMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1),
      mutatorPtr_(0),
      outputFileName_(),
      initialWeightFileName_("0"),
      stepCount_(0)
   {  setClassName("WangLandauAdaptiveStepMove"); } 
   
   /* 
   * Read parameter speciesId.
   */
   void WangLandauAdaptiveStepMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<GeneralpolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Error: Species must be GeneralpolymerSG");
      }
      read<Pair <int> >(in, "Range", Range_);
      if (Range_[1]-Range_[0] < 0) {
         UTIL_THROW("Error: Total range is negative");
      }
      Species* speciesPtr = &system().simulation().species(speciesId_);
      capacity_ = speciesPtr->capacity()+1;
       
      mutatorPtr_ = &speciesPtr_->mutator();
      weights_.allocate(capacity_);
      stateCount_.allocate(capacity_);      
      read<double>(in, "weightStep", weightSize_);
      read<std::string>(in, "outputFileName",outputFileName_);
      readOptional<std::string>(in, "initialWeights",initialWeightFileName_);
      std::ifstream weightFile;
      if (initialWeightFileName_!="0") {
         system().fileMaster().open(initialWeightFileName_, weightFile);
         int n;
         double m;
         while (weightFile >> n >>m)
         {
           weights_[n]= m;
         }
      } else {
         for (int x = 0; x < capacity_; ++x) {
             weights_[x]=0;
         }
        }
      std::string fileName = outputFileName_;
      for (int x = 0; x < capacity_; ++x) {
        stateCount_[x] = 0; 
      }
        fileName = outputFileName_+".weights";      
        system().fileMaster().openOutputFile(fileName, outputFile_);
        outputFile_ << stepCount_ << "	" << weightSize_ << std::endl;
        crit_=.1;
   }
   /*
   * Load state from an archive.
   */
   void WangLandauAdaptiveStepMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);

      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<GeneralpolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Species is not a GeneralpolymerSG");
      }
      loadParameter<Pair <int> >(ar, "Range", Range_);
      if (Range_[1]-Range_[0] < 0) {
         UTIL_THROW("Error: Total range is negative");
      }
      weights_.allocate(capacity_);      
      loadParameter<double>(ar, "weightStep", weightSize_); 
   }
   

   /*
   * Save state to an archive.
   */
   void WangLandauAdaptiveStepMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;  
   }
   // Determine if it is time to adapt the step size and if so do such
   void WangLandauAdaptiveStepMove::stepAdapt()
   { 
     // Determine if the histogram is sufficiently flat
     bool flat = true;
     int movesMade = 0;
     for (int y = Range_[0]; y <= Range_[1]; ++y) {
       movesMade = movesMade + stateCount_[y];
     }
     double binAve = movesMade/(Range_[1]-Range_[0]+1);
     
     for (int z = Range_[0]; z <= Range_[1]; ++z) {
       if (std::abs(stateCount_[z]-binAve)/binAve > crit_) {
       flat = false;
       }
     }

     // If so Adapt and clear the histogram
     if (flat) {
        crit_=crit_/2;
        weightSize_ = pow(weightSize_, .5);
        outputFile_ << stepCount_ << "	    " << weightSize_ << std::endl;
       for (int x = 0; x < capacity_; ++x) {
           stateCount_[x] = 0;
       }
     }

   }

   Molecule& WangLandauAdaptiveStepMove::randomSGMolecule(int speciesId, int nType0, int flipType)
   {
      int moleculeId,nType,nMol,index,type;
      int count = 0;
      nMol = system().nMolecule(speciesId);
      if (nMol <= 0) {
         Log::file() << "Number of molecules in species " << speciesId
                     << " = " << nMol << std::endl;
         UTIL_THROW("Number of molecules in species <= 0");
      }
   /*   if (nType0==Range_[0]) {
         flipType = 1;
         nType=nMol-nType0;
      } else {
      if (nType0==Range_[1]) {
         flipType = 0;
         nType=nType0;
      } else {
         return system().randomMolecule(speciesId);
      }
      } 
     */ 
      if (flipType == 0) {
         nType = nType0;
      } else 
      if (flipType == 1) {
         nType = nMol-nType0;
      } else 
      if (flipType ==2) {
      
      return system().randomMolecule(speciesId);
      }
      index = simulation().random().uniformInt(0, nType);
      for (int i=0; i<nMol; ++i) {
        type = speciesPtr_->mutator().moleculeStateId(system().molecule(speciesId, i));
         if (type==flipType) {
            if (count==index) {
               moleculeId = i;
               return system().molecule(speciesId, moleculeId);
            }
            count = count + 1;
         }
      }
      UTIL_THROW("No molecule selected");
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool WangLandauAdaptiveStepMove::move() 
   {  stepCount_=stepCount_+1;
      incrementNAttempt();
      // Special semigrand selector
      int oldState = mutatorPtr_->stateOccupancy(0);
      int nMolecule = system().nMolecule(speciesId_);
      if (oldState == Range_[1] || oldState == Range_[0]) {
         if (oldState == Range_[1]) {
            if (simulation().random().uniformInt(1,nMolecule+1) > oldState)
            {
            bool accept = false;
            int state = mutatorPtr_->stateOccupancy(0);
            weights_[state]=weights_[state]+log(weightSize_);
            stateCount_[state]=stateCount_[state]+1;
            stepAdapt();
            return accept;
       
            }
            else {
              flipType_=0;}
         }
         if (oldState == Range_[0]) {
            if (simulation().random().uniformInt(1,nMolecule+1) <= oldState) {
            bool accept =false;
            int state = mutatorPtr_->stateOccupancy(0);
            weights_[state]=weights_[state]+log(weightSize_);
            stateCount_[state]=stateCount_[state]+1;
            incrementNAccept();
            return accept;
               
    } else {
              flipType_=1;}  
        }
      } else
      if (oldState > Range_[1]) {
         flipType_=0;
         Molecule& molecule = randomSGMolecule(speciesId_, oldState,flipType_);
         speciesPtr_->mutator().setMoleculeState(molecule,1);
         bool accept = true;
         incrementNAccept();
         return accept;
      } else
      if (oldState < Range_[0]) {
         flipType_=1;
         Molecule& molecule = randomSGMolecule(speciesId_, oldState,flipType_);
         speciesPtr_->mutator().setMoleculeState(molecule,0);
         bool accept = true;
         incrementNAccept();
         return accept;
      } else {
        flipType_=2;
      }

      Molecule& molecule = randomSGMolecule(speciesId_, oldState, flipType_);
      #ifndef SIMP_NOPAIR
      // Calculate pair energy for the chosen molecule
      double oldEnergy = system().pairPotential().moleculeEnergy(molecule);
      #endif
      // Toggle state of the molecule
      int oldStateId = speciesPtr_->mutator().moleculeStateId(molecule);
      int newStateId = (oldStateId == 0) ? 1 : 0;
      speciesPtr_->mutator().setMoleculeState(molecule, newStateId);
      int stateChange = -1;
      if  (newStateId == 0) {
          stateChange = 1;
      }

      #ifdef SIMP_NOPAIR 

      bool   accept = true;
      bool   inAllowedRange = true;
      
      #else // ifndef SIMP_NOPAIR

      // Recalculate pair energy for the molecule
      double newEnergy = system().pairPotential().moleculeEnergy(molecule);

      // Decide whether to accept or reject
      int    newState = mutatorPtr_->stateOccupancy(0);
      // Different move if the move is with in the desired range or not
      //int    oldState = newState - stateChange;
      double oldWeight = weights_[oldState];
      double newWeight = weights_[newState];
      double ratio  = boltzmann(newEnergy - oldEnergy)*exp(oldWeight-newWeight);
      //double ratio = exp(oldWeight-newWeight);
      bool   accept = random().metropolis(ratio);
      #endif

      if (accept) {

         incrementNAccept();
      } else {

         // Revert chosen molecule to original state
         speciesPtr_->mutator().setMoleculeState(molecule, oldStateId);
      }
      int state = mutatorPtr_->stateOccupancy(0);
      weights_[state]=weights_[state]+log(weightSize_);
      stateCount_[state]=stateCount_[state]+1;
      stepAdapt();
      return accept;
   }
 
   void WangLandauAdaptiveStepMove::output()
   {    outputFile_.close();
       
        std::string fileName = outputFileName_; 
        std::ofstream outputFile;
        fileName += ".dat";
        system().fileMaster().openOutputFile(fileName, outputFile);
        for (int i = 0; i < capacity_; i++) {
           outputFile << i << "   " <<  weights_[i]<<std::endl;
        }
        outputFile.close();
   }
   

}
