/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WangLandauAdaptiveStepMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/species/GeneralpolymerSG.h>
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
      mutatorPtr_ = &speciesPtr_->mutator();
      weights_.allocate(Range_[1]-Range_[0]+1);
      stateCount_.allocate(Range_[1]-Range_[0]+1);      
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
      }   
      for (int x = 0; x < Range_[1]-Range_[0]+1; ++x) {
          stateCount_[x] = 0;
      }
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
      weights_.allocate(Range_[1]-Range_[0]+1);      
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
     for (int y = 0; y < Range_[1]-Range_[0]+1; ++y) {
       movesMade = movesMade + stateCount_[y];
     }
     double binAve = movesMade/(Range_[1]-Range_[0]+1);
     
     for (int z = 0; z < Range_[1] - Range_[0]+1; ++z) {
       if (std::abs(stateCount_[z]-binAve)/binAve > 0.07) {
       flat = false;
       }
     }

     // If so Adapt and clear the histogram
     if (flat) {
        weightSize_ = pow(weightSize_, .5);
        weightTrack_[i]=weightSize_;
        steps_[i]=stepCount_    
       for (int x = 0; x < Range_[1]-Range_[0]+1; ++x) {
           stateCount_[x] = 0;
       }
     }

   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool WangLandauAdaptiveStepMove::move() 
   {  stepCount_=stepCount_+1
      incrementNAttempt();
      Molecule& molecule = system().randomMolecule(speciesId_);

      #ifndef INTER_NOPAIR
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

      #ifdef INTER_NOPAIR 

      bool   accept = true;

      #else // ifndef INTER_NOPAIR

      // Recalculate pair energy for the molecule
      double newEnergy = system().pairPotential().moleculeEnergy(molecule);

      // Decide whether to accept or reject
      int    newState = mutatorPtr_->stateOccupancy(0);
      int    oldState = newState - stateChange;
      double oldWeight = weights_[oldState-Range_[0]];
      double newWeight = weights_[newState-Range_[0]];
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
      weights_[state-Range_[0]]=weights_[state-Range_[0]]+log(weightSize_);
      stateCount_[state-Range_[0]]=stateCount_[state-Range_[0]]+1;
      stepAdapt();
      return accept;
   }
 
   void WangLandauAdaptiveStepMove::output()
   {
        std::ofstream outputFile;
        std::string fileName = outputFileName_;
        fileName += ".dat";
        system().fileMaster().openOutputFile(fileName, outputFile);
        for (int i = 0; i < Range_[1]-Range_[0]+1; i++) {

           outputFile << i+Range_[0] << "   " <<  weights_[i]<<std::endl;
        }
        outputFile.close();
        // File of the time steps and weights
        fileName += ".steps";
        system().fileMaster().openOutputFile(fileName, outputFile);
        for (int i = 0; i < weightTrack_.size(); i++) {

           outputFile << steps_[i] << "		" <<  weightTrack_[i]<<std::endl;
        }
        outputFile.close();
   }
   

}
