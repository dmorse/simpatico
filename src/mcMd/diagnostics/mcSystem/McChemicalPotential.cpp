#ifndef MCMD_MC_CHEMICAL_POTENTIAL_CPP
#define MCMD_MC_CHEMICAL_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McChemicalPotential.h"                        // class header
#include <util/misc/FileMaster.h>  
#include <util/archives/Serializable_includes.h>

#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcMoves/SystemMove.h>
#include <mcMd/mcSimulation/mc_potentials.h>

#include <util/boundary/Boundary.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McChemicalPotential::McChemicalPotential(McSystem& system)
    : SystemDiagnostic<McSystem>(system),
      systemPtr_(&system),
      simulationPtr_(&system.simulation()),
      boundaryPtr_(&system.boundary()),
      randomPtr_(&system.simulation().random()),
      outputFile_(),
      accumulator_(),
      nTrial_(-1),
      nMoleculeTrial_(-1),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {
      energyEnsemblePtr_ = &system.energyEnsemble();
   }

   /*
   * Read parameters and initialize.
   */
   void McChemicalPotential::readParam(std::istream& in)
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      read<int>(in, "nTrial", nTrial_);
      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }

      read<int>(in, "nMoleculeTrial", nMoleculeTrial_);
      if (nMoleculeTrial_ <=0 || nMoleculeTrial_ > MaxMoleculeTrial_) {
         UTIL_THROW("Invalid value input for nMoleculeTrial");
      }

      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ <=0 || speciesId_ > system().simulation().nSpecies()) {
         UTIL_THROW("Invalid value input for nMoleculeTrial");
      }

      isInitialized_ = true;
   }

   /*
   * Clear accumulator.
   */
   void McChemicalPotential::setup() 
   {  accumulator_.clear(); }
 
   
   /* 
   * Evaluate Rosenbluth weight, and add to accumulator.
   */
   void McChemicalPotential::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         accumulator_.sample(chemicalPotential(), outputFile_);
      }     
   }


   /*
   * Output results to file after simulation is completed.
   */
   void McChemicalPotential::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();
   }
   
   
   /*
   * Save state to binary file archive.
   */
   void McChemicalPotential::save(Serializable::OArchive& ar)
   {  ar & *this; }
   
   /*
   * Load state from a binary file archive.
   */
   void McChemicalPotential::load(Serializable::IArchive& ar)
   { ar & *this; }
   

   /*
   * Configuration bias algorithm for adding one atom to a chain end.
   */
   void
   McChemicalPotential::addEndAtom(Atom* endPtr, Atom* pvtPtr, int bondType,
                          double &rosenbluth, double &energy)
   {
      Vector trialPos[MaxTrial_];
      Vector bondVec;
      Vector pvtPos = pvtPtr->position();
      double trialProb[MaxTrial_], trialEnergy[MaxTrial_];
      double beta, length;
      int    iTrial;

      #ifdef INTER_ANGLE
      Species::AtomAngleArray angles;
      const Angle *anglePtr;
      const Atom  *pvtPtr2(NULL);
      Vector dr1, dr2;
      int    iAngle, angleTypeId(0);
      double rsq1, rsq2, cosTheta;
      #endif
   
      // Generate a random bond length
      beta   = energyEnsemble().beta();
      length = 
         system().bondPotential().randomBondLength(&random(), beta, bondType);
   
      // Loop over nTrial trial positions:
      rosenbluth = 0.0;
      for (iTrial=0; iTrial < nTrial_; ++iTrial) {
         random().unitVector(bondVec);
         bondVec *= length;
         // trialPos = pvtPos + bondVec
         trialPos[iTrial].add(pvtPos, bondVec); 
         boundary().shift(trialPos[iTrial]);
         endPtr->position() = trialPos[iTrial];
         #ifndef INTER_NOPAIR
         trialEnergy[iTrial] = system().pairPotential().atomEnergy(*endPtr);
         #else
         trialEnergy[iTrial] = 0.0;
         #endif

         #ifdef INTER_ANGLE
         if (system().hasAnglePotential()) {

            endPtr->molecule().species().getAtomAngles(*endPtr, angles);
            for (iAngle = 0; iAngle < angles.size(); ++iAngle) {
               anglePtr = angles[iAngle];
               if (&anglePtr->atom(1) == pvtPtr) {
                  if (&anglePtr->atom(0) == endPtr) {
                     pvtPtr2 = &anglePtr->atom(2);
                  } else {
                     pvtPtr2 = &anglePtr->atom(0);
                  }
                  angleTypeId = anglePtr->typeId();
               }
            }
   
            // Get the angle spanned.
            rsq1 = boundary().distanceSq(pvtPtr->position(),
                                         pvtPtr2->position(), dr1);
            rsq2 = boundary().distanceSq(endPtr->position(),
                                         pvtPtr->position(), dr2);
            cosTheta = dr1.dot(dr2) / sqrt(rsq1 * rsq2);
   
            trialEnergy[iTrial] += system().anglePotential().
                                   energy(cosTheta, angleTypeId);
         }
         #endif

         #ifdef INTER_EXTERNAL
         trialEnergy[iTrial] += system().externalPotential().atomEnergy(*endPtr);
         #endif

         trialProb[iTrial] = boltzmann(trialEnergy[iTrial]);
         rosenbluth += trialProb[iTrial];
      }
   
      // Normalize trial probabilities 
      for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
         trialProb[iTrial] = trialProb[iTrial]/rosenbluth;
      }
     
      // Choose trial position
      iTrial = random().drawFrom(trialProb, nTrial_);
   
      // Calculate total energy for chosen trial.
      energy = system().bondPotential().energy(length*length, bondType);
      energy += trialEnergy[iTrial];
   
      // Set position of new end atom to chosen value
      endPtr->position() = trialPos[iTrial];

   }


   /*
   * Configuration bias algorithm for deleting one atom from chain end.
   */
   void 
   McChemicalPotential::deleteEndAtom(Atom* endPtr, Atom* pvtPtr, int bondType,
                            double &rosenbluth, double &energy)
   {
      Vector bondVec;
      Vector pvtPos = pvtPtr->position();
      double trialEnergy, lengthSq, length;
   
      // Calculate bond length of pvt-end bond
      lengthSq = boundary().distanceSq(pvtPos, endPtr->position());
      length   = sqrt(lengthSq);

      // Calculate current nonbonded pair energy of end atom
      #ifndef INTER_NOPAIR
      energy = system().pairPotential().atomEnergy(*endPtr);
      #else
      energy = 0.0;
      #endif

      #ifdef INTER_ANGLE
      Species::AtomAngleArray angles;
      const Angle *anglePtr;
      const Atom  *pvtPtr2(NULL);
      Vector dr1, dr2;
      int    iAngle, angleTypeId(0);
      double rsq1, rsq2, cosTheta;

      if (system().hasAnglePotential()) {

         // Get the angle type and pointers of atoms forming the angle.
         endPtr->molecule().species().getAtomAngles(*endPtr, angles);
         for (iAngle = 0; iAngle < angles.size(); ++iAngle) {
            anglePtr = angles[iAngle];
            if (&anglePtr->atom(1) == pvtPtr) {
               if (&anglePtr->atom(0) == endPtr) {
                  pvtPtr2 = &anglePtr->atom(2);
               } else {
                  pvtPtr2 = &anglePtr->atom(0);
               }
               angleTypeId = anglePtr->typeId();
            }
         }
   
         // Calculate angle energy. 
         rsq1 = boundary().distanceSq(pvtPtr->position(),
                                      pvtPtr2->position(), dr1);
         rsq2 = boundary().distanceSq(endPtr->position(),
                                      pvtPtr->position(), dr2);
         cosTheta = dr1.dot(dr2) / sqrt(rsq1 * rsq2);
   
         energy += system().anglePotential().energy(cosTheta, angleTypeId);

      }
      #endif

      #ifdef INTER_EXTERNAL
      if (system().hasExternalPotential()) {
         energy += system().externalPotential().atomEnergy(*endPtr);
      }
      #endif

      // Rosenbluth factor = exp(-beta*(pair + angle + external))
      rosenbluth = boltzmann(energy);
      
      // Add bond energy of current pvt-end bond to energy.
      // This is the final value of inout energy parameter.
      energy += system().bondPotential().energy(lengthSq, bondType);

      // Loop over nTrial - 1 additional trial positions:
      for (int iTrial=0; iTrial < nTrial_ - 1; ++iTrial) {

         random().unitVector(bondVec);
         bondVec *= length;
         endPtr->position().add(pvtPos, bondVec);  
         boundary().shift(endPtr->position());

         #ifndef INTER_NOPAIR
         trialEnergy = system().pairPotential().atomEnergy(*endPtr);
         #else
         trialEnergy = 0.0;
         #endif

         #ifdef INTER_ANGLE
         if (system().hasAnglePotential()) {

            // Get the angle type and atom pointer at the angle.
            endPtr->molecule().species().getAtomAngles(*endPtr, angles);
            for (iAngle = 0; iAngle < angles.size(); ++iAngle) {
               anglePtr = angles[iAngle];
               if (&anglePtr->atom(1) == pvtPtr) {
                  if (&anglePtr->atom(0) == endPtr) {
                     pvtPtr2 = &anglePtr->atom(2);
                  } else {
                     pvtPtr2 = &anglePtr->atom(0);
                  }
                  angleTypeId = anglePtr->typeId();
               }
            }

            // Calculate energy.
            rsq1 = boundary().distanceSq(pvtPtr->position(),
                                         pvtPtr2->position(), dr1);
            rsq2 = boundary().distanceSq(endPtr->position(),
                                         pvtPtr->position(), dr2);
            cosTheta = dr1.dot(dr2) / sqrt(rsq1 * rsq2);
   
            trialEnergy += system().anglePotential()
                                   .energy(cosTheta, angleTypeId);
         }
         #endif

         #ifdef INTER_EXTERNAL
         if (system().hasExternalPotential()) {
            trialEnergy += system().externalPotential()
                                   .atomEnergy(*endPtr);
         }
         #endif

         rosenbluth += boltzmann(trialEnergy);
      }

   }

   /*
   * Rosenbluth algorithm for calculationg chemical potential.
   */
   double McChemicalPotential::chemicalPotential() 
   { 
      Species* speciesPtr;
      Molecule* molPtr;
      Molecule::BondIterator bondIter;
      Atom* endPtr;
      double rosenbluth;
      double e = 0;
 
      int nSpecies; 
      nSpecies = system().simulation().nSpecies();
         if (nSpecies != 1) {
            UTIL_THROW("Error: nSpecies != 1");
         }

      speciesPtr = &(simulation().species(speciesId_));
      molPtr = &(speciesPtr->reservoir().pop());
      system().addMolecule(*molPtr);

      endPtr = &molPtr->atom(0);
      boundary().randomPosition(random(), endPtr->position());
      
      system().pairPotential().addAtom(*endPtr);
      rosenbluth = boltzmann(system().pairPotential().atomEnergy(*endPtr));

      for (molPtr->begin(bondIter); bondIter.notEnd(); ++bondIter) {
          addEndAtom(&(bondIter->atom(1)), &(bondIter->atom(0)), bondIter->typeId(), rosenbluth, e);
          rosenbluth *= rosenbluth;
          e += e;
          system().pairPotential().addAtom(bondIter->atom(1));
      }

      rosenbluth = rosenbluth/pow(nTrial_,molPtr->nAtom());

      for (molPtr->begin(bondIter); bondIter.notEnd(); ++bondIter) {
          system().pairPotential().deleteAtom(bondIter->atom(1));
      }
      
      system().removeMolecule(*molPtr);

      return rosenbluth; 
   }

}
#endif 
