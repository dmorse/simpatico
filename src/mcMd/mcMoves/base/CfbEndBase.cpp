/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbEndBase.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcSimulation/mc_potentials.h>

#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/getAtomGroups.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   CfbEndBase::CfbEndBase(McSystem& system) : 
      SystemMove(system),
      nTrial_(-1)
   {
      // Precondition
      #ifdef SIMP_DIHEDRAL
      if (system.hasDihedralPotential()) {
         UTIL_THROW("CfbEndBase is unusable with dihedrals");
      }
      #endif
   }
  
   /* 
   * Destructor
   */
   CfbEndBase::~CfbEndBase() 
   {}
   
   /* 
   * Read parameter nTrial.
   */
   void CfbEndBase::readParameters(std::istream& in) 
   {
      read<int>(in, "nTrial", nTrial_);
      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }
   }
 
   /*
   * Configuration bias algorithm for deleting one atom from chain end.
   */
   void 
   CfbEndBase::deleteEndAtom(Atom* endPtr, Atom* pvtPtr, int bondType,
                            double &rosenbluth, double &energy)
   {
      Vector bondVec;
      Vector pvtPos = pvtPtr->position();
      double trialEnergy, lengthSq, length;
   
      // Calculate bond length of pvt-end bond
      lengthSq = boundary().distanceSq(pvtPos, endPtr->position());
      length   = sqrt(lengthSq);

      // Calculate current nonbonded pair energy of end atom
      #ifndef SIMP_NOPAIR
      energy = system().pairPotential().atomEnergy(*endPtr);
      #else
      energy = 0.0;
      #endif

      #ifdef SIMP_ANGLE
      AtomAngleArray angles;
      const Angle *anglePtr;
      const Atom  *pvtPtr2(NULL);
      Vector dr1, dr2;
      int    iAngle, angleTypeId(0);
      double rsq1, rsq2, cosTheta;

      if (system().hasAnglePotential()) {

         // Get the angle type and pointers of atoms forming the angle.
         getAtomAngles(*endPtr, angles);
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

      #ifdef SIMP_EXTERNAL
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

         #ifndef SIMP_NOPAIR
         trialEnergy = system().pairPotential().atomEnergy(*endPtr);
         #else
         trialEnergy = 0.0;
         #endif

         #ifdef SIMP_ANGLE
         if (system().hasAnglePotential()) {

            // Get the angle type and atom pointer at the angle.
            getAtomAngles(*endPtr, angles);
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

         #ifdef SIMP_EXTERNAL
         if (system().hasExternalPotential()) {
            trialEnergy += 
               system().externalPotential().atomEnergy(*endPtr);
         }
         #endif

         rosenbluth += boltzmann(trialEnergy);
      }

   }
   
  
   /*
   * Configuration bias algorithm for adding one atom to a chain end.
   */
   void
   CfbEndBase::addEndAtom(Atom* endPtr, Atom* pvtPtr, int bondType,
                          double &rosenbluth, double &energy)
   {
      Vector trialPos[MaxTrial_];
      Vector bondVec;
      Vector pvtPos = pvtPtr->position();
      double trialProb[MaxTrial_], trialEnergy[MaxTrial_];
      double beta, length;
      int    iTrial;

      #ifdef SIMP_ANGLE
      AtomAngleArray angles;
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
         #ifndef SIMP_NOPAIR
         trialEnergy[iTrial] = system().pairPotential().atomEnergy(*endPtr);
         #else
         trialEnergy[iTrial] = 0.0;
         #endif

         #ifdef SIMP_ANGLE
         if (system().hasAnglePotential()) {

            getAtomAngles(*endPtr, angles);
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

         #ifdef SIMP_EXTERNAL
         if (system().hasExternalPotential()) {
            trialEnergy[iTrial] += 
                system().externalPotential().atomEnergy(*endPtr);
         }
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

}
