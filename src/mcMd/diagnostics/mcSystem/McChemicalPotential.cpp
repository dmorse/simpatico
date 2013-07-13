#ifndef MCMD_MC_CHEMICAL_POTENTIAL_CPP
#define MCMD_MC_CHEMICAL_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McChemicalPotential.h"                        // class header

#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcMoves/SystemMove.h>
#include <mcMd/mcSimulation/mc_potentials.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>

#include <util/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
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
      energyEnsemblePtr_(&system.energyEnsemble()),
      randomPtr_(&system.simulation().random()),
      outputFile_(),
      accumulator_(),
      nTrial_(-1),
      nMoleculeTrial_(-1),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {}

   /*
   * Read parameters and initialize.
   */
   void McChemicalPotential::readParameters(std::istream& in)
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      read<int>(in, "nTrial", nTrial_);
      read<int>(in, "nMoleculeTrial", nMoleculeTrial_);
      read<int>(in, "speciesId", speciesId_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }

      if (nMoleculeTrial_ <=0) {
         UTIL_THROW("Invalide value input for nMoleculeTrial");
      }

      if (speciesId_ <0 || speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("Invalid value input for speciesId");
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

         Species* speciesPtr;
         Molecule* molPtr;
         Molecule::BondIterator bondIter;
         Atom* endPtr;
         double w;
         double rosenbluth = 1;
         double de;
         double e = 0;

         speciesPtr = &(simulation().species(speciesId_));

         // Pop a new molecule off the species reservoir
         molPtr = &(speciesPtr->reservoir().pop());
         system().addMolecule(*molPtr);

         // Loop over molecule growth trials
         for (int i = 0; i < nMoleculeTrial_; i++) {

            // Pick a random position for the first atom
            endPtr = &molPtr->atom(0);
            boundary().randomPosition(random(), endPtr->position());

            e = system().pairPotential().atomEnergy(*endPtr);
            rosenbluth = boltzmann(e);
            system().pairPotential().addAtom(*endPtr);

            for (molPtr->begin(bondIter); bondIter.notEnd(); ++bondIter) {
                addEndAtom(&(bondIter->atom(1)), &(bondIter->atom(0)), bondIter->typeId(), w, de);
                e += de;
                rosenbluth *= w;
                system().pairPotential().addAtom(bondIter->atom(1));
            }

            rosenbluth = rosenbluth / pow(nTrial_,molPtr->nAtom()-1);
            accumulator_.sample(rosenbluth, outputFile_);

            system().pairPotential().deleteAtom(*endPtr);
            for (molPtr->begin(bondIter); bondIter.notEnd(); ++bondIter) {
                system().pairPotential().deleteAtom(bondIter->atom(1));
            }
         }

         // Return additional molecule to reservoir
         system().removeMolecule(*molPtr);
         speciesPtr->reservoir().push(*molPtr);
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

      // Open and write a *.prm file for the parameter block
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      // Open a *.ave file for Average accumulator output
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

      // Calculate total energy for chosen trial
      energy = system().bondPotential().energy(length*length, bondType);
      energy += trialEnergy[iTrial];

      // Set position of new end atom to chosen trialPos Vector
      endPtr->position() = trialPos[iTrial];
   }

}
#endif
