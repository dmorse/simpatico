/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WangLandauMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#ifndef SIMP_NOPAIR
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
   WangLandauMove::WangLandauMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1),
      mutatorPtr_(0),
      outputFileName_()
   {  setClassName("WangLandauMove"); } 
   
   /* 
   * Read parameter speciesId.
   */
   void WangLandauMove::readParameters(std::istream& in) 
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
      read<double>(in, "weightStep", weightSize_);
      for (int x = 0; x < Range_[1]-Range_[0]+1; ++x) {
          weights_[x] = 0;
      }
      read<std::string>(in, "outputFileName",outputFileName_);
   }
   /*
   * Load state from an archive.
   */
   void WangLandauMove::loadParameters(Serializable::IArchive& ar)
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
   void WangLandauMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;  
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool WangLandauMove::move() 
   {
      incrementNAttempt();
      Molecule& molecule = system().randomMolecule(speciesId_);

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

      #else // ifndef SIMP_NOPAIR

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
      return accept;
   }
 
   void WangLandauMove::output()
   {
        std::ofstream outputFile;
        std::string fileName = outputFileName_;
        fileName += ".dat";
        system().fileMaster().openOutputFile(fileName, outputFile);
        for (int i = 0; i < Range_[1]-Range_[0]+1; i++) {

           outputFile << i+Range_[0] << "   " <<  weights_[i]<<std::endl;
        }
        outputFile.close();
   }
   

}
