/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McNVTChemicalPotential.h"                        // class header
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>

#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcMoves/SystemMove.h>
#include <mcMd/mcSimulation/mc_potentials.h>

#include <util/boundary/Boundary.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/global.h>

#include <math.h>
#include <cstdio>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McNVTChemicalPotential::McNVTChemicalPotential(McSystem& system)
    : SystemAnalyzer<McSystem>(system),
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
   void McNVTChemicalPotential::readParameters(std::istream& in)
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      read<int>(in, "nTrial", nTrial_);
      read<int>(in, "nMoleculeTrial", nMoleculeTrial_);
      read<int>(in, "speciesId", speciesId_);
      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      read<double>(in, "Emin", Emin_);
      read<double>(in, "Emax", Emax_);
      read<double>(in, "EnBin", EnBin_);
      Eaccumulator_.setParam(Emin_, Emax_, EnBin_);
      Eaccumulator_.clear();
      read<double>(in, "Ecmin", Ecmin_);
      read<double>(in, "Ecmax", Ecmax_);
      read<double>(in, "EcnBin", EcnBin_);
      Ecaccumulator_.setParam(Ecmin_, Ecmax_, EcnBin_);
      Ecaccumulator_.clear();
      read<double>(in, "Emmin", Emmin_);
      read<double>(in, "Emmax", Emmax_);
      read<double>(in, "EmnBin", EmnBin_);
      Emaccumulator_.setParam(Emmin_, Emmax_, EmnBin_);
      Emaccumulator_.clear();
      read<double>(in, "BRmin", BRmin_);
      read<double>(in, "BRmax", BRmax_);
      read<double>(in, "BRnBin", BRnBin_);
      BRaccumulator_.setParam(BRmin_, BRmax_, BRnBin_);
      BRaccumulator_.clear();

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }

      if (nMoleculeTrial_ <=0) {
         UTIL_THROW("Invalid value input for nMoleculeTrial");
      }

      if (speciesId_ <0 || speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("Invalid value input for speciesId");
      }

      isInitialized_ = true;

   }

   /*
   * Clear accumulator.
   */
   void McNVTChemicalPotential::setup()
   {  accumulator_.clear(); }

   /*
   * Evaluate Rosenbluth weight, and add to accumulator.
   */
   void McNVTChemicalPotential::sample(long iStep)
   {
      if (isAtInterval(iStep))  {

         Molecule* molPtr;
         Atom* endPtr, *pvt1Ptr;
         Vector trialPos[MaxTrial_], bondVec, pvt1Pos;
         #ifdef SIMP_ANGLE
         double angle;
         Atom* pvt2Ptr;
         Vector pvt2Pos;
         #endif
         double beta, length;
         // Rosenbluth factor and energy of inserted bonds at each stage.
         double w = 0, de = 0;
         // Total rosenbluth factor and energy of inserted polymer.
         double rosenbluth = 1, energy = 0;
         double trialProb[MaxTrial_], trialEnergy[MaxTrial_];
         int    iTrial;

         beta = energyEnsemble().beta();
         molPtr = &(simulation().getMolecule(speciesId_));
         system().addMolecule(*molPtr);

         // Main molecule inserting loop:
         for (int i = 0; i < nMoleculeTrial_; i++) {

            // Inserting the first atom.
            endPtr = &molPtr->atom(0);
            for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
               // trialPos = pvtPos + bondVec
               boundary().randomPosition(random(), trialPos[iTrial]);
               boundary().shift(trialPos[iTrial]);

               #ifndef SIMP_NOPAIR
               trialEnergy[iTrial] = 
                  system().pairPotential().atomEnergy(*endPtr);
               #else
               trialEnergy[iTrial] = 0.0;
               #endif

               #ifdef SIMP_EXTERNAL
               trialEnergy[iTrial] += 
                  system().externalPotential().atomEnergy(*endPtr);
               #endif

               trialProb[iTrial] = boltzmann(trialEnergy[iTrial]);
               w += trialProb[iTrial];
            }

            // Normalize trial probabilities
            for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
               trialProb[iTrial] = trialProb[iTrial]/rosenbluth;
            }

            // Choose trial position
            iTrial = random().drawFrom(trialProb, nTrial_);

            // Calculate total energy for chosen trial
            de = trialEnergy[iTrial];

            // Set position of new end atom to chosen trialPos Vector
            endPtr->position() = trialPos[iTrial];

            #ifndef SIMP_NOPAIR
            system().pairPotential().addAtom(*endPtr);
            #endif

            rosenbluth *= w;
            energy += de;

            // Inserting the second atom.
            int bondTypeId = molPtr->bond(0).typeId();
            length = system().bondPotential().randomBondLength(&random(), 
                                                               beta, 
                                                               bondTypeId);

            w = 0;
            de = 0;

            pvt1Ptr = &molPtr->atom(0);
            endPtr = &molPtr->atom(1);

            pvt1Pos = pvt1Ptr->position();

            for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
               // trialPos = pvtPos + bondVec
               random().unitVector(bondVec);
               bondVec *= length;
               // trialPos = pvtPos + bondVec
               trialPos[iTrial].add(pvt1Pos, bondVec);
               boundary().shift(trialPos[iTrial]);
               endPtr->position() = trialPos[iTrial];

               #ifndef SIMP_NOPAIR
               trialEnergy[iTrial] =
                  system().pairPotential().atomEnergy(*endPtr);
               #else
               trialEnergy[iTrial] = 0.0;
               #endif

               #ifdef SIMP_EXTERNAL
               trialEnergy[iTrial] += 
                  system().externalPotential().atomEnergy(*endPtr);
               #endif

               trialProb[iTrial] = boltzmann(trialEnergy[iTrial]);
               w += trialProb[iTrial];
            }

            // Normalize trial probabilities
            for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
               trialProb[iTrial] = trialProb[iTrial]/rosenbluth;
            }

            // Choose trial position
            iTrial = random().drawFrom(trialProb, nTrial_);

            // Fetch nonbonded energy for chosen trial
            de = trialEnergy[iTrial];

            // Add bond energy
            de += system().bondPotential().energy(length*length, bondTypeId);

            // Set position of new end atom to chosen trialPos Vector
            endPtr->position() = trialPos[iTrial];

            #ifndef SIMP_NOPAIR
            system().pairPotential().addAtom(*endPtr);
            #endif

            rosenbluth *= w;
            energy += de;

            // Inserting the third atom.
            bondTypeId = molPtr->bond(1).typeId();
            length = 
               system().bondPotential().randomBondLength(&random(), 
                                                         beta, bondTypeId);

            #ifdef SIMP_ANGLE
            int angleTypeId = molPtr->angle(0).typeId();
            angle = 
               system().anglePotential().randomAngle(&random(), 
                                                     beta, angleTypeId);
            pvt2Ptr = &molPtr->atom(0);
            pvt2Pos = pvt2Ptr->position();
            Vector n = pvt1Pos;
            n -= pvt2Pos;
            #endif

            w = 0;
            de = 0;
            pvt1Ptr = &molPtr->atom(1);
            pvt1Pos = pvt1Ptr->position();
            endPtr = &molPtr->atom(2);

            for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
               // trialPos = pvtPos + bondVec
               #ifdef SIMP_ANGLE
               if (system().hasAnglePotential()) {
                  uniformCone(length, angle, n, bondVec);
               } else {
                  random().unitVector(bondVec);
                  bondVec *= length;
               }
               #else
               random().unitVector(bondVec);
               bondVec *= length;
               #endif

               // trialPos = pvtPos + bondVec
               trialPos[iTrial].add(pvt1Pos, bondVec);
               boundary().shift(trialPos[iTrial]);
               endPtr->position() = trialPos[iTrial];

               #ifndef SIMP_NOPAIR
               trialEnergy[iTrial] = 
                  system().pairPotential().atomEnergy(*endPtr);
               #else
               trialEnergy[iTrial] = 0.0;
               #endif

               #ifdef SIMP_EXTERNAL
               trialEnergy[iTrial] += 
                  system().externalPotential().atomEnergy(*endPtr);
               #endif

               trialProb[iTrial] = boltzmann(trialEnergy[iTrial]);
               w += trialProb[iTrial];
            }
            // Normalize trial probabilities
            for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
               trialProb[iTrial] = trialProb[iTrial]/rosenbluth;
            }

            // Choose trial position
            iTrial = random().drawFrom(trialProb, nTrial_);

            // Get nonbonded energy for chosen trial
            de = trialEnergy[iTrial];

            // Add bonded energies
            de += system().bondPotential().energy(length*length, bondTypeId);
            #ifdef SIMP_ANGLE
            de += system().anglePotential().energy(cos(angle), angleTypeId);
            #endif

            // Set position of new end atom to chosen trialPos Vector
            endPtr->position() = trialPos[iTrial];

            #ifndef SIMP_NOPAIR
            system().pairPotential().addAtom(*endPtr);
            #endif

            rosenbluth *= w;
            energy += de;

            for (int atomId = 3; atomId < molPtr->nAtom(); ++atomId) {
                addEndAtom(molPtr, atomId, w, de);
                rosenbluth *= w;
                energy += de;
                #ifndef SIMP_NOPAIR
                system().pairPotential().addAtom(molPtr->atom(atomId));
                #endif
            }

            rosenbluth = rosenbluth / pow(nTrial_,molPtr->nAtom());
            accumulator_.sample(rosenbluth, outputFile_);

            #ifndef SIMP_NOPAIR
            for (int atomId = 0; atomId < molPtr->nAtom(); ++atomId) {
                system().pairPotential().deleteAtom(molPtr->atom(atomId));
            }
            #endif
         }

         system().removeMolecule(*molPtr);
         simulation().returnMolecule(*molPtr);

         }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void McNVTChemicalPotential::output()
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

      fileMaster().openOutputFile(outputFileName(".e"), outputFile_);
      Eaccumulator_.output(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".ec"), outputFile_);
      Ecaccumulator_.output(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".emin"), outputFile_);
      Emaccumulator_.output(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".br"), outputFile_);
      BRaccumulator_.output(outputFile_);
      outputFile_.close();
   }

   /*
   * Save state to binary file archive.
   */
   void McNVTChemicalPotential::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void McNVTChemicalPotential::load(Serializable::IArchive& ar)
   { ar & *this; }


    /*
   * Configuration bias algorithm for adding one atom to a chain end.
   */
   void
   McNVTChemicalPotential::addEndAtom(Molecule* molPtr, 
                                      int atomId, double &rosenbluth, 
                                      double &energy)
   {
      Atom* endPtr = &molPtr->atom(atomId);
      Vector trialPos[MaxTrial_];
      Vector bondVec;
      Vector pvt1Pos = molPtr->atom(atomId-1).position();

      #ifdef SIMP_ANGLE
      Vector pvt2Pos = molPtr->atom(atomId-2).position();
      #endif

      #ifdef SIMP_DIHEDRAL
      Vector pvt3Pos = molPtr->atom(atomId-3).position();
      #endif

      double trialProb[MaxTrial_], trialEnergy[MaxTrial_];
      double beta, length;
      int bondTypeId = molPtr->bond(atomId - 1).typeId();

      #ifdef SIMP_ANGLE
      double angle;
      int angleTypeId = molPtr->angle(atomId - 2).typeId();
      #endif

      int    iTrial;
      double eMin = 0;

      // Generate a random bond-length and bond-angle.
      beta   = energyEnsemble().beta();
      length = system().bondPotential().randomBondLength(&random(), beta, 
                                        molPtr->bond(atomId-1).typeId());
      #ifdef SIMP_ANGLE
      Vector n = pvt1Pos;
      n -= pvt2Pos;
      angle = system().anglePotential().randomAngle(&random(), beta, 
                                        molPtr->angle(atomId-2).typeId());
      #endif
      #ifdef SIMP_DIHEDRAL
      Vector dR1 = pvt3Pos;
      Vector dR2 = pvt2Pos;
      Vector dR3 = pvt1Pos;
      dR1 -= pvt2Pos;
      dR2 -= pvt1Pos;
      dR3 -= endPtr->position();
      #endif

      // Loop over nTrial trial positions:
      rosenbluth = 0.0;
      for (iTrial=0; iTrial < nTrial_; ++iTrial) {

         #ifdef SIMP_ANGLE
         if (system().hasAnglePotential()) {
            uniformCone(length, angle, n, bondVec);
         } else {
            random().unitVector(bondVec);
            bondVec *= length;
         }
         #else
         random().unitVector(bondVec);
         bondVec *= length;
         #endif

         trialPos[iTrial].add(pvt1Pos, bondVec);
         boundary().shift(trialPos[iTrial]);
         endPtr->position() = trialPos[iTrial];

         #ifndef SIMP_NOPAIR
         trialEnergy[iTrial] = system().pairPotential().atomEnergy(*endPtr);
         #else
         trialEnergy[iTrial] = 0.0;
         #endif

         #ifdef SIMP_DIHEDRAL
         if (system().hasDihedralPotential()) {
            trialEnergy[iTrial] += 
               system().dihedralPotential().energy(dR1, dR2, dR3, 
                                        molPtr->dihedral(atomId-3).typeId());
         }
         #endif

         #ifdef SIMP_EXTERNAL
         trialEnergy[iTrial] += 
            system().externalPotential().atomEnergy(*endPtr);
         #endif

         // Finding the Minimum of the trialEnergy vector
         if (eMin > trialEnergy[iTrial])
            eMin = trialEnergy[iTrial];
         Eaccumulator_.sample(trialEnergy[iTrial]);

         trialProb[iTrial] = boltzmann(trialEnergy[iTrial]);
         rosenbluth += trialProb[iTrial];
      }

      Emaccumulator_.sample(eMin);
      BRaccumulator_.sample(rosenbluth);
      // Normalize trial probabilities
      for (iTrial = 0; iTrial < nTrial_; ++iTrial) {
         trialProb[iTrial] = trialProb[iTrial]/rosenbluth;
      }

      // Choose trial position
      iTrial = random().drawFrom(trialProb, nTrial_);

      // Calculate total energy for chosen trial.
      energy = system().bondPotential().energy(length*length, bondTypeId);
      #ifdef SIMP_ANGLE
      energy += system().anglePotential().energy(cos(angle), angleTypeId);
      #endif

      energy += trialEnergy[iTrial];
      Ecaccumulator_.sample(energy);

      // Set position of new end atom to chosen value
      endPtr->position() = trialPos[iTrial];
   }

   /*
   * Finding a vector which make an angle \theta with vector n and has length b.
   */
   void McNVTChemicalPotential::uniformCone(const double length, 
                                            const double angle, 
                                            const Vector n, Vector &p)
   {
      Vector v1;

      do {
            random().unitVector(p);
            v1.divide(n,n.abs());
            v1.multiply(v1,p.projection(v1));
            p -= v1;
      } while (p.abs() > 1.0E-8);

      v1.divide(n,n.abs());
      v1 *= length*cos(angle);
      p /= p.abs();
      p *= length*sin(angle);
      p += (v1);
   }

}
