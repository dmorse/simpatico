/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WangLandauALTMove.h"
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
   WangLandauALTMove::WangLandauALTMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1),
      mutatorPtr_(0),
      outputFileName_(),
      initialWeightFileName_("0"),
      stepCount_(0)
   {  setClassName("WangLandauALTMove"); } 
   
   /* 
   * Read parameter speciesId.
   */
   void WangLandauALTMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<GeneralpolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Error: Species must be GeneralpolymerSG");
      }
      Species* speciesPtr = &system().simulation().species(speciesId_);
      capacity_ = speciesPtr->capacity()+1; 
      mutatorPtr_ = &speciesPtr_->mutator();
      weights_.allocate(capacity_);
      stateCount_.allocate(capacity_); 
      read<double>(in, "weightStep", weightSize_);
      read<int>(in, "UpperLimit", Ulimit_);
      read<int>(in, "LowerLimit", Llimit_);
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
        crit_=.05;
   }
   /*
   * Load state from an archive.
   */
   void WangLandauALTMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);

      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<GeneralpolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Species is not a GeneralpolymerSG");
      }
      weights_.allocate(capacity_);      
      loadParameter<double>(ar, "weightStep", weightSize_); 
      loadParameter<int>(ar, "UpperLimit", Ulimit_);
      loadParameter<int>(ar, "LowerLimit", Llimit_);
   }
   

   /*
   * Save state to an archive.
   */
   void WangLandauALTMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;  
   }
   // Determine if it is time to adapt the step size and if so do such
   void WangLandauALTMove::stepAdapt()
   { 
     // Determine if the histogram is sufficiently flat
     bool flat = true;
     int movesMade = 0;
     for (int y = Llimit_; y <= Ulimit_; ++y) {
       movesMade = movesMade + stateCount_[y];
     }
     double binAve = movesMade/(Ulimit_-Llimit_+1);
     
     for (int z = Llimit_; z <= Ulimit_; ++z) {
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

   Molecule& WangLandauALTMove::randomSGMolecule(int speciesId, int nSubType, int flipType)
   {
      int moleculeId,nMol,index,type;
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
       
      index = simulation().random().uniformInt(0, nSubType);
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
   bool WangLandauALTMove::move() 
   {  stepCount_=stepCount_+1;
      incrementNAttempt();
      // Special semigrand selector
           SpeciesMutator* mutatorPtr = &speciesPtr_->mutator();
      int oldStateCount = mutatorPtr->stateOccupancy(0);
      if  (oldStateCount == Ulimit_) {
      flipSubtype_ = simulation().random().uniformInt(0,2);
      //flipSubtype_ = 0;
      if (flipSubtype_ == 1) {
        bool accept = false;
        return accept;
        }
      }
      if (oldStateCount == Llimit_) {
      flipSubtype_ = simulation().random().uniformInt(0,2);
      //flipSubtype_=1;
      if (flipSubtype_ == 0){
         bool accept = false;
         return accept;
       }
      }
      if (oldStateCount>Llimit_ && oldStateCount < Ulimit_) {
      flipSubtype_ = simulation().random().uniformInt(0,2);
      }
      int oldState = mutatorPtr->stateOccupancy(0);
      oldStateCount = mutatorPtr->stateOccupancy(flipSubtype_);
      Molecule& molecule = randomSGMolecule(speciesId_, oldStateCount, flipSubtype_);
      #ifndef SIMP_NOPAIR
      // Calculate pair energy for the chosen molecule
      double oldEnergy = system().pairPotential().moleculeEnergy(molecule);
      #endif
      // Toggle state of the molecule
      int oldStateId = speciesPtr_->mutator().moleculeStateId(molecule);
      int newStateId = (oldStateId == 0) ? 1 : 0;
      speciesPtr_->mutator().setMoleculeState(molecule, newStateId);
  //    int stateChange = -1;
   //   if  (newStateId == 0) {
     //     stateChange = 1;
     // }

      #ifdef SIMP_NOPAIR 

      bool   accept = true;
    //  bool   inAllowedRange = true;
      
       #else //ifndef SIMP_NOPAIR

      // Recalculate pair energy for the molecule
      double newEnergy = system().pairPotential().moleculeEnergy(molecule);
      double numoWeight = (double)(mutatorPtr->stateOccupancy(oldStateId)+1)/(double)(mutatorPtr->stateOccupancy(newStateId));
      // Decide whether to accept or reject
      int    newState = mutatorPtr_->stateOccupancy(0);
      // Different move if the move is with in the desired range or not
      //int    oldState = newState - stateChange;
      double oldWeight = weights_[oldState];
      double newWeight = weights_[newState];
      double ratio  = boltzmann(newEnergy - oldEnergy)*exp(oldWeight-newWeight)*numoWeight;
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
 
   void WangLandauALTMove::output()
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
