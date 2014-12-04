/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbLinear.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcSimulation/mc_potentials.h>

#include <mcMd/chemistry/getAtomGroups.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>

#include <util/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   CfbLinear::CfbLinear(McSystem& system) :
      SystemMove(system),
      nTrial_(-1)
   {
      // Precondition
      #ifdef INTER_DIHEDRAL
      if (system.hasDihedralPotential()) {
         UTIL_THROW("CfbLinear is unusable with dihedrals");
      }
      #endif
   }

   /*
   * Destructor
   */
   CfbLinear::~CfbLinear()
   {}

   /*
   * Read parameter nTrial.
   */
   void CfbLinear::readParameters(std::istream& in)
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
   CfbLinear::deleteAtom(Molecule& molecule, Atom& atom0, int atomId, int sign,
                            double &rosenbluth, double &energy)
   {
      // sign == 0, direction == -1  -> atomId = 0, 1,... end
      // sign == 1, direction == +1 -> atomId = nAtom-1, nAtom-2, etc.
      assert(sign == 0 || sign == 1);
      int direction = sign ? 1 : -1;
      Vector& pos0 = atom0.position();
      Vector& pos1 = molecule.atom(atomId - direction).position();

      // Calculate bond length and energy
      Vector v1, u1;
      double r1, bondEnergy;
      int bondType;
      r1 = boundary().distanceSq(pos1, pos0, v1); // v1 = pos1 - pos0
      bondType = molecule.bond(atomId - sign).typeId();
      bondEnergy = system().bondPotential().energy(r1, bondType);
      r1 = sqrt(r1); // bond length
      u1 = v1;
      u1 /= r1;      // unit vector

      // Calculate other initial energies (excluding bond)

      #ifndef INTER_NOPAIR
      energy = system().pairPotential().atomEnergy(atom0);
      #else
      energy = 0.0;
      #endif

      #ifdef INTER_ANGLE
      Vector v2, u2;
      double r2, cosTheta;
      Vector* pos2Ptr;
      int angleType;
      if (molecule.nAngle()) {
         pos2Ptr = molecule.atom(atomId - 2*direction).position();
         r2 = boundary().distanceSq(*pos2Ptr, pos1, v2); // v2 = pos2 - pos1
         r2 = sqrt(r2);
         u2 = v2;
         u2 /= r2; // unit vector
         double cosTheta = u1.dot(u2);
         angleType = molecule.bond(atomId - 2*sign).typeId();
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

      // Add bond energy to total energy in current position
      energy += bondEnergy;

      // Loop over nTrial - 1 additional trial positions:
      double trialEnergy;
      for (int iTrial = 0; iTrial < nTrial_ - 1; ++iTrial) {

         // Generate trial position
         random().unitVector(u1);
         v1 = u1;
         v1 *= r1;
         pos0.subtract(pos1, v1);
         boundary().shift(pos0);

         // Evaluate trial energy (excluding bond energy)
         #ifndef INTER_NOPAIR
         trialEnergy = system().pairPotential().atomEnergy(atom0);
         #else
         trialEnergy = 0.0;
         #endif
         #ifdef INTER_ANGLE
         if (molecule.nAngle()) {
            cosTheta = u1.dot(u2);
            trialEnergy += system().anglePotential()
                                   .energy(cosTheta, angleType);
         }
         #endif
         #ifdef INTER_EXTERNAL
         if (system().hasExternalPotential()) {
            trialEnergy += system().externalPotential().atomEnergy(atom0);
         }
         #endif

         rosenbluth += boltzmann(trialEnergy);
      }

   }

   /*
   * Configuration bias algorithm for adding one atom to a chain end.
   */
   void
   CfbLinear::addAtom(Molecule& molecule, Atom& atom0, int atomId, int sign,
                      double &rosenbluth, double &energy)
   {
      Vector trialPos[MaxTrial_], v1, u1;
      double trialProb[MaxTrial_], trialEnergy[MaxTrial_];

      int direction = sign ? 1 : -1;
      Vector& pos0 = atom0.position();
      Vector& pos1 = molecule.atom(atomId - direction).position();

      // Generate a random bond length
      double beta, r1;
      int bondType, iTrial;
      beta = energyEnsemble().beta();
      r1 = system().bondPotential().randomBondLength(&random(), beta, bondType);
      bondType = molecule.bond(atomId - sign).typeId();

      #ifdef INTER_ANGLE
      // Calculate vector v2 = pos2 - pos1, r2 = |v2|
      Vector v2, u2;
      double r2, cosTheta;
      Vector* pos2Ptr;
      int angleType;
      if (molecule().nAngle()) {
         angleType = molecule.bond(atomId - 2*sign).typeId();
         pos2Ptr = molecule.atom(atomId - 2*direction).position();
         r2 = boundary().distanceSq(*pos2Ptr, pos1, v2);
         r2 = sqrt(r2);
         u2 = v2;
         u2 /= r2;
      }
      #endif

      // Loop over nTrial trial positions:
      rosenbluth = 0.0;
      for (iTrial = 0; iTrial < nTrial_; ++iTrial) {

         // Generate trial position
         random().unitVector(u1);
         v1 = u1;
         v1 *= r1;
         trialPos[iTrial].subtract(pos1, v1);
         boundary().shift(trialPos[iTrial]);
         pos0 = trialPos[iTrial];

         // Compute trial energy (excluding bond energy)
         #ifndef INTER_NOPAIR
         trialEnergy[iTrial] = system().pairPotential().atomEnergy(atom0);
         #else
         trialEnergy[iTrial] = 0.0;
         #endif
         #ifdef INTER_ANGLE
         if (molecule().nAngle()) {
            cosTheta = u1.dot(u2);
            trialEnergy[iTrial] += system().anglePotential().
                                   energy(cosTheta, angleTypeId);
         }
         #endif
         #ifdef INTER_EXTERNAL
         if (system().hasExternalPotential()) {
            trialEnergy[iTrial] +=
                        system().externalPotential().atomEnergy(atom0);
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
      energy = trialEnergy[iTrial];
      energy += system().bondPotential().energy(r1*r1, bondType);

      // Set position of end atom to chosen value
      pos0 = trialPos[iTrial];
   }

}
