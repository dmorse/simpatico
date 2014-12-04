/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
      #ifdef INTER_DIHEDRAL
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
   CfbEndBase::deleteAtom(Molecule& molecule, int atomId, int sign,
                            double &rosenbluth, double &energy)
   {
      int direction = sign ? 1 : -1;
      Atom& atom0 = molecule.atom(atomId);
      Vector& pos0 = atom0.position();

      // Calculate length of end bond
      Atom& atom1 = molecule.atom(atomId - direction);
      Vector& pos1 = atom1.position();
      Vector v1;
      double rsq1 = boundary().distanceSq(pos1, pos0, v1);
      double r1 = sqrt(rsq1);
      int bondType = molecule.bond(atomId - sign).typeId();

      // Calculate current nonbonded pair energy of end atom
      #ifndef INTER_NOPAIR
      energy = system().pairPotential().atomEnergy(atom0);
      #else
      energy = 0.0;
      #endif

      #ifdef INTER_ANGLE
      Vector v2;
      Vector* pos2Ptr;
      int angleType;
      double rsq2, r2, cosTheta;

      if (molecule.nAngle()) {
         pos2Ptr = molecule.atom(atomId - 2*direction).position();
         rsq2 = boundary().distanceSq(*pos2Ptr, *pos1Ptr, v2);
         r2 = sqrt(rsq2);
         cosTheta = v1.dot(v2) / (r1*r2);
         int angleType = molecule.bond(atomId - 2*sign).typeId();
         energy += system().anglePotential().energy(cosTheta, angleType);
      }
      #endif

      #ifdef INTER_EXTERNAL
      if (system().hasExternalPotential()) {
         energy += system().externalPotential().atomEnergy(atom0);
      }
      #endif

      // Rosenbluth factor = exp(-beta*(pair + angle + external))
      rosenbluth = boltzmann(energy);

      // Add current bond energy to energy.
      energy += system().bondPotential().energy(rsq1, bondType);

      // Loop over nTrial - 1 additional trial positions:
      Vector bondVec;
      for (int iTrial=0; iTrial < nTrial_ - 1; ++iTrial) {

         random().unitVector(v1);
         v1 *= r1;
         pos0.subtract(pos1, v1);
         boundary().shift(pos0);

         #ifndef INTER_NOPAIR
         trialEnergy = system().pairPotential().atomEnergy(atom0);
         #else
         trialEnergy = 0.0;
         #endif

         #ifdef INTER_ANGLE
         if (molecule.nAngle()) {

            cosTheta = v1.dot(v2) / (r1*r2);
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
   * Configuration bias algorithm for adding one atom to a chain end.
   */
   void
   CfbEndBase::addAtom(Molecule& molecule, int atomId, int sign,
                       double &rosenbluth, double &energy)
   {
      int direction = sign ? 1 : -1;
      Atom& atom0 = molecule.atom(atomId);
      Atom& atom1 = molecule.atom(atomId - direction);
      Vector& pos0 = atom0.position();
      Vector& pos1 = atom1.position();

      // Generate a random bond length
      beta = energyEnsemble().beta();
      r1 = system().bondPotential().randomBondLength(&random(), beta, bondType);
      Vector trialPos[MaxTrial_];
      Vector v1;
      double trialProb[MaxTrial_], trialEnergy[MaxTrial_];
      double beta, length;
      int iTrial;

      #ifdef INTER_ANGLE
      Vector v2;
      int angleType;
      double rsq2, r2, cosTheta;
      #endif

      // Loop over nTrial trial positions:
      rosenbluth = 0.0;
      for (iTrial=0; iTrial < nTrial_; ++iTrial) {
         random().unitVector(v1);
         bondVec *= length;
         trialPos[iTrial].subtract(pos1, v1);
         boundary().shift(trialPos[iTrial]);
         pos0 = trialPos[iTrial];
         #ifndef INTER_NOPAIR
         trialEnergy[iTrial] = system().pairPotential().atomEnergy(atom0);
         #else
         trialEnergy[iTrial] = 0.0;
         #endif

         #ifdef INTER_ANGLE
         if (system().hasAnglePotential()) {
            cosTheta = v1.dot(v2) / (r1*r2);
            trialEnergy[iTrial] += system().anglePotential().
                                   energy(cosTheta, angleTypeId);
         }
         #endif

         #ifdef INTER_EXTERNAL
         trialEnergy[iTrial] +=
                        system().externalPotential().atomEnergy(atom0);
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
      energy = system().bondPotential().energy(r1*r1, bondType);
      energy += trialEnergy[iTrial];

      // Set position of new end atom to chosen value
      pos0 = trialPos[iTrial];

   }

}
