/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GeneralpolymerALTSemiGrandMove.h"
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
   GeneralpolymerALTSemiGrandMove::GeneralpolymerALTSemiGrandMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1)
   {  setClassName("GeneralpolymerALTSemiGrandMove"); } 
   
   /* 
   * Read parameter speciesId.
   */
   void GeneralpolymerALTSemiGrandMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "UpperLimit", Ulimit_);
      read<int>(in, "LowerLimit", Llimit_);
      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<GeneralpolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Error: Species must be GeneralpolymerSG");
      }
  
   }
 
   /*
   * Load state from an archive.
   */
   void GeneralpolymerALTSemiGrandMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);

      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<GeneralpolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Species is not a GeneralpolymerSG");
      }
   }

   /*
   * Save state to an archive.
   */
   void GeneralpolymerALTSemiGrandMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;  
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool GeneralpolymerALTSemiGrandMove::move() 
   {
      incrementNAttempt();
      // choose whether to flip type 0 or type 1
      SpeciesMutator* mutatorPtr = &speciesPtr_->mutator();
      int oldStateCount = mutatorPtr->stateOccupancy(0);
      if  (oldStateCount == Ulimit_) {
      flipSubtype_ = simulation().random().uniformInt(0,2);
      if (flipSubtype_ == 1) {
        bool accept = false;
        return accept;
        }
      } 
      if (oldStateCount == Llimit_) {
      flipSubtype_ = simulation().random().uniformInt(0,2);
      if (flipSubtype_ == 0){
         bool accept = false;
         return accept;
       }
      } 
      if (oldStateCount>Llimit_ && oldStateCount < Ulimit_) {
      flipSubtype_ = simulation().random().uniformInt(0,2);
      }
       // Instead of flipping a random molecule, flip one using the selective flipper
      //Molecule& molecule = system().randomMolecule(speciesId_);
      oldStateCount = mutatorPtr->stateOccupancy(flipSubtype_);
      Molecule& molecule = randomSGMolecule(speciesId_, oldStateCount, flipSubtype_);
      //SpeciesMutator* mutatorPtr = &speciesPtr_->mutator();
      #ifndef SIMP_NOPAIR
      // Calculate pair energy for the chosen molecule
      double oldEnergy = system().pairPotential().moleculeEnergy(molecule);
      #endif

      // Toggle state of the molecule
      int oldStateId = speciesPtr_->mutator().moleculeStateId(molecule);
      int newStateId = (oldStateId == 0) ? 1 : 0;
      speciesPtr_->mutator().setMoleculeState(molecule, newStateId);
      #ifdef SIMP_NOPAIR 
      bool   accept = true;
      #else // ifndef SIMP_NOPAIR
      int newStateTotal = mutatorPtr->stateOccupancy(0);
      // Recalculate pair energy for the molecule
      double newEnergy = system().pairPotential().moleculeEnergy(molecule);
      double numoWeight = (double)(mutatorPtr->stateOccupancy(oldStateId)+1)/(double)(mutatorPtr->stateOccupancy(newStateId));
      // Decide whether to accept or reject
      double oldWeight = speciesPtr_->mutator().stateWeight(oldStateId);
      double newWeight = speciesPtr_->mutator().stateWeight(newStateId);
      double ratio  = boltzmann(newEnergy - oldEnergy)*newWeight/oldWeight*numoWeight;
      bool   accept = random().metropolis(ratio);
      #endif

      if (accept) {
      
           incrementNAccept();

      } else {

         // Revert chosen molecule to original state
         speciesPtr_->mutator().setMoleculeState(molecule, oldStateId);
         bool accept = false;
         return accept;
      }

      return accept;
   }
 Molecule& GeneralpolymerALTSemiGrandMove::randomSGMolecule(int speciesId, int nSubType,  int flipType)
   {
      int moleculeId,nMol,index,type;
      int count = 0;
      //int shouldFlip = 1;
      nMol = system().nMolecule(speciesId);
      if (nMol <= 0) {
         Log::file() << "Number of molecules in species " << speciesId
                     << " = " << nMol << std::endl;
         UTIL_THROW("Number of molecules in species <= 0");
      }
      index = simulation().random().uniformInt(0, nSubType);
      for (int i=0; i<nMol; ++i) {
        type = speciesPtr_->mutator().moleculeStateId(system().molecule(speciesId, i));
//        std::cout << type << "     " << moleculeId << "     ";
         if (type==flipType) {
            if (count==index) {
               moleculeId = i;
               //std::cout << type << "     " << moleculeId << "     "; 
               return system().molecule(speciesId, moleculeId);
            }
          count = count+1;
         }
      }
     UTIL_THROW("Nothing Selected");
   }


}
