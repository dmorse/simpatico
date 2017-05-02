/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbRebridgeBase.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcSimulation/mc_potentials.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/chemistry/Bond.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   CfbRebridgeBase::CfbRebridgeBase(McSystem& system) : 
      CfbEndBase(system),
      length21_(1.0),
      length10_(1.0),
      kappa10_(0.0)
   {
      // Preconditions
      #ifdef MCMD_NOMASKBONDED
      UTIL_THROW("CfbRebridgeBase is unusable ifdef MCMD_NOMASKBONDED");
      #endif
      #ifdef SIMP_DIHEDRAL
      if (system.hasDihedralPotential()) {
         UTIL_THROW("CfbRebrigeBase is unusable with dihedrals");
      }
      #endif
   } 
   /* 
   * Destructor.
   */
   CfbRebridgeBase::~CfbRebridgeBase() 
   {} 
   
   /* 
   * Read parameters from file.
   */
   void CfbRebridgeBase::readParameters(std::istream& in) 
   {
      // Number of trial configuration bias attempts.
      read<int>(in, "nTrial", nTrial_);
      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }

      // Parameters needed by orientation bias.
      read<double>(in, "length21", length21_);
      read<double>(in, "length10", length10_);
      read<double>(in, "kappa10", kappa10_);
   }

   /* 
   * Load state from archive.
   */
   void CfbRebridgeBase::loadParameters(Serializable::IArchive& ar) 
   {
      loadParameter<int>(ar, "nTrial", nTrial_);
      loadParameter<double>(ar, "length21", length21_);
      loadParameter<double>(ar, "length10", length10_);
      loadParameter<double>(ar, "kappa10", kappa10_);

      // Validate
      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }
   }

   /* 
   * Save state to archive.
   */
   void CfbRebridgeBase::save(Serializable::OArchive& ar) 
   {
      ar & nTrial_;
      ar & length21_;
      ar & length10_;
      ar & kappa10_;
   }

   /* 
   * Initialize tables for spring constants and normalizations.
   */
   void CfbRebridgeBase::setup()
   {
      double l20;
      double lmin, lmax, del, pi;
      int    i;
   
      pi = acos(-1.0);
   
      // Range of l_20; 2.0 is empirical.
      lmax = 2.0 * (length10_ + length21_);
   
      // Fill table elements for: l_20 <= lmin
      if (length21_ == length10_) {
         // 0.1 is empirical.
         lmin = length10_ / length21_ / sqrt(kappa10_) * 0.1;
         l20Table[0]    = lmin;
         angTable[0]    = pi/2.0;
         kappaTable[0]  = 0.0;
         normalTable[0] = 2.0;
      } else {
         // 0.5 is empirical.
         lmin = fabs(length10_ - length21_) * 0.5;
         l20Table[0] = lmin;
         orientationBiasTable(lmin, angTable[0], kappaTable[0], normalTable[0]);
      }
   
      // Fill table elements for: lmin < l_20 < lmax
      del  = (lmax - lmin) / double(MaxBin_-2);
      l20  = lmin - del/2.0;
      for (i=1; i<MaxBin_-1; i++) {
         l20 += del;
         l20Table[i] = l20;
         orientationBiasTable(l20, angTable[i], kappaTable[i], normalTable[i]);
      }
   
      // Fill table elements for: l_20 >= lmax
      l20Table[MaxBin_-1]    = lmax;
      angTable[MaxBin_-1]    = 0.0;
      kappaTable[MaxBin_-1]  = kappa10_ * length21_ * lmax;
      computeNormalization(angTable[MaxBin_-1], kappaTable[MaxBin_-1], 
                           normalTable[MaxBin_-1]);
   }

   /* 
   * Configuration bias algorithm for deleting the last particle of an interior
   * bridge. 
   *
   *                     1 
   *                      o-----o 0
   *                     /
   *                    o
   *                    2
   *         2 = prev, 0 = next, 1 = crank-shaft rotating
   *
   * Note that in the current implementation, the end atom position will be
   * modified, so the calling subroutine is responsible for recovering the old
   * position if necessary. This is similar to the implementation of
   * CfbEndBase::deleteEndAtom.
   */
   void CfbRebridgeBase::deleteMiddleAtom(
                Atom* partPtr, Atom* prevPtr, Atom* nextPtr,
                int prevBType, int nextBType,
                double &rosenbluth, double &energy)
   {
      Vector  prevPos = prevPtr->position();
      Vector  nextPos = nextPtr->position();
      Vector  partPos = partPtr->position();
      Vector  bondVec, trialPos, u_20, u_21;
      double  trialEnergy, bondEnergy, wExt, length, lengthSq;
      double  l_20, l_21, bias;
      double  prefAng, kappaAng, normConst;
      int     iTrial;
      bool    ready;
      
      // Bond length and bonding vector 2 --> 0
      lengthSq = boundary().distanceSq(nextPos, prevPos, u_20);
      l_20  = sqrt(lengthSq);
      u_20 /= l_20;
   
      // Bond length and bonding vector 2 --> 1 & bonding energy
      lengthSq = boundary().distanceSq(partPos, prevPos, u_21);
      l_21  = sqrt(lengthSq);
      u_21 /= l_21;
      energy = system().bondPotential().energy(lengthSq, prevBType);
   
      // Bond length between 1 and 0 & bonding energy
      lengthSq   = boundary().distanceSq(partPos, nextPos);
      bondEnergy = system().bondPotential().energy(lengthSq, nextBType);
      energy    += bondEnergy;
   
      // Nonbonded pair energy of current position 
      #ifndef SIMP_NOPAIR
      trialEnergy = system().pairPotential().atomEnergy(*partPtr);
      #else
      trialEnergy = 0.0;
      #endif

      #ifdef SIMP_ANGLE
      if (system().hasAnglePotential()) {
         trialEnergy += system().anglePotential().atomEnergy(*partPtr);
      }
      #endif

      #ifdef SIMP_EXTERNAL
      if (system().hasExternalPotential()) {
         trialEnergy += system().externalPotential().atomEnergy(*partPtr);
      }
      #endif

      energy += trialEnergy;
    
      // Compute orientation/Angle bias factor and the normalization factor.
      getAngKappaNorm(l_20, prefAng, kappaAng, normConst);
      orientationBias(u_21, u_20, prefAng, kappaAng, bias);
   
      // remove normalized orientation bias
      rosenbluth = normConst / bias;
   
      // start to accumulate Rosenbluth factor
      wExt = boltzmann(trialEnergy + bondEnergy);
      
      // Generate "nTrial-1" trial positions based on 1-0 bonding.
      length = l_21;
      for (iTrial=0; iTrial < nTrial_ - 1; iTrial++) {
         ready = false;
         while (!ready) {
            random().unitVector(bondVec);
            bondVec *= length;
            trialPos.add(prevPos, bondVec);
            boundary().shift(trialPos);
   
            lengthSq = boundary().distanceSq(trialPos, prevPos, u_21);
            l_21  = sqrt(lengthSq);
            u_21 /= l_21;
            orientationBias(u_21, u_20, prefAng, kappaAng, bias);
   
   	    // Rejection method
   	    if (bias > random().uniform(0.0, 1.0)) {
               ready = true;
            }
         }
         
         // Bond 1-0 potential energy
         lengthSq   = boundary().distanceSq(trialPos, nextPos);
         bondEnergy = system().bondPotential().energy(lengthSq, nextBType);
         partPtr->position() = trialPos;
         #ifndef SIMP_NOPAIR
         trialEnergy = system().pairPotential().atomEnergy(*partPtr);
         #else
         trialEnergy = 0.0;
         #endif

         #ifdef SIMP_ANGLE
         if (system().hasAnglePotential()) {
            trialEnergy += system().anglePotential().atomEnergy(*partPtr);
         }
         #endif

         #ifdef SIMP_EXTERNAL
         if (system().hasExternalPotential()) {
            trialEnergy += system().externalPotential().atomEnergy(*partPtr);
         }
         #endif

         wExt += boltzmann(trialEnergy + bondEnergy);
      }
   
      // Update Rosenbluth weight (orientation bias has been removed earlier)
      rosenbluth *= wExt;
   
   }
   
   /* 
   * Configuration bias algorithm for adding the last particle of an interior
   * bridge. 
   */
   void CfbRebridgeBase::addMiddleAtom(
                 Atom* partPtr, Atom* prevPtr, Atom* nextPtr,
                 int prevBType, int nextBType,
                 double &rosenbluth, double &energy)
   {
      Vector  prevPos = prevPtr->position();
      Vector  nextPos = nextPtr->position();
      Vector  u_20, u_21, bondVec;
      Vector  trialPos[MaxTrial_];
      double  trialEnergy[MaxTrial_], bondEnergy[MaxTrial_];
      double  trialProb[MaxTrial_], bias[MaxTrial_];
      double  l_20, l_21;
      double  beta, length, lengthSq, wExt;
      double  prefAng, kappaAng, normConst;
      int     iTrial;
      bool    ready;

      // Bond 2 --> 0
      lengthSq = boundary().distanceSq(nextPos, prevPos, u_20);
      l_20     = sqrt(lengthSq);
      u_20 /= l_20;
 
      // Bond 2-1 & bonding energy
      beta   = energyEnsemble().beta();
      length = system().bondPotential().
               randomBondLength(&random(), beta, prevBType);
      l_21   = length;
      energy = system().bondPotential().energy(l_21*l_21, prevBType);
   
      // Loop over trial orientations
      wExt = 0.0;
      
      // Compute normalization factor for the orientation bias.
      getAngKappaNorm(l_20, prefAng, kappaAng, normConst);
   
      // Generating trial orientations
      for (iTrial=0; iTrial < nTrial_; iTrial++) {
         ready = false;
         while (!ready) {
            random().unitVector(bondVec);
            bondVec *= length;
            trialPos[iTrial].add(prevPos, bondVec);
            boundary().shift(trialPos[iTrial]);
   
            lengthSq = boundary().distanceSq(trialPos[iTrial], prevPos, u_21);
            l_21  = sqrt(lengthSq);
            u_21 /= l_21;
            orientationBias(u_21, u_20, prefAng, kappaAng, bias[iTrial]);
   
   	    // Rejection method
   	    if (bias[iTrial] > random().uniform(0.0, 1.0)) {
               ready = true;
            }
         }
         
         // Bond 1-0 potential energy 
         lengthSq = boundary().distanceSq(trialPos[iTrial], nextPos);
         bondEnergy[iTrial] = 
            system().bondPotential().energy(lengthSq, nextBType);
         partPtr->position() = trialPos[iTrial];

         #ifndef SIMP_NOPAIR
         trialEnergy[iTrial] = system().pairPotential().atomEnergy(*partPtr);
         #else
         trialEnergy[iTrial] = 0.0;
         #endif

         #ifdef SIMP_ANGLE
         if (system().hasAnglePotential()) {
            trialEnergy[iTrial] += 
                  system().anglePotential().atomEnergy(*partPtr);
         }
         #endif

         #ifdef SIMP_EXTERNAL
         if (system().hasExternalPotential()) {
            trialEnergy[iTrial] += 
                  system().externalPotential().atomEnergy(*partPtr);
         }
         #endif

         trialProb[iTrial] =
            boltzmann(trialEnergy[iTrial] + bondEnergy[iTrial]);
         wExt += trialProb[iTrial];
      }
   
      // Normalize trial probabilities 
      for (iTrial=0; iTrial < nTrial_; iTrial++) {
         trialProb[iTrial] = trialProb[iTrial] / wExt;
      }
     
      // Choose trial position
      iTrial = random().drawFrom(trialProb, nTrial_);
   
      // Add non-bonded pair energy for chosen trial
      energy += trialEnergy[iTrial];
      
      // Add 1-0 bond energy
      energy += bondEnergy[iTrial];
      
      // Copy trial position to newPos
      partPtr->position() = trialPos[iTrial];
   
      // update Rosenbluth weight and remove the orientation bias
      rosenbluth = wExt * normConst / bias[iTrial];
   }


   /* 
   * Delete a consequtive sequence of atoms.
   */
   void CfbRebridgeBase::deleteSequence(int nActive, int sign, Atom* endPtr,
              int *bonds, double &rosenbluth, double &energy)
   {
      Atom   *thisPtr;
      Atom   *prevPtr;
      Atom   *nextPtr;
      int     prevBType, nextBType;
      double  rosen_r, energy_r;

      // Initialize weight
      rosenbluth = 1.0;
      energy = 0.0;

      thisPtr = endPtr;
      prevPtr = thisPtr - sign;
      nextPtr = thisPtr + sign;
      // Step 1: delete the last atom
      if (nActive >= 1) {
         nextBType = bonds[0];
         prevBType = bonds[1];
 
         // Orientation biased trimer rebridge move
         deleteMiddleAtom(thisPtr, prevPtr, nextPtr,
                prevBType, nextBType, rosen_r, energy_r);
         #ifndef SIMP_NOPAIR
         system().pairPotential().deleteAtom(*thisPtr);
         #endif
         rosenbluth *= rosen_r;
         energy += energy_r;
      }

      // step 2: delete the remaining atoms
      for (int i = 0; i < nActive - 1; ++i) {
         thisPtr -= sign;
         prevPtr -= sign;

         prevBType = bonds[i+2];
         deleteEndAtom(thisPtr, prevPtr, prevBType, rosen_r, energy_r);
         #ifndef SIMP_NOPAIR
         system().pairPotential().deleteAtom(*thisPtr);
         #endif

         rosenbluth *= rosen_r;
         energy += energy_r;
      }

   }

   /* 
   * Add a consequtive sequence of atoms.
   */
   void CfbRebridgeBase::addSequence(int nActive, int sign, Atom* beginPtr,
              int *bonds, double &rosenbluth, double &energy)
   {
      Atom   *thisPtr;
      Atom   *prevPtr;
      Atom   *nextPtr;
      int     prevBType, nextBType;
      double  rosen_f, energy_f;

      // Initialize weight
      rosenbluth = 1.0;
      energy = 0.0;

      // Step 1: grow the first nActive - 1 atoms
      thisPtr = beginPtr;
      prevPtr = thisPtr - sign;

      for (int i = 0; i < nActive - 1; ++i) {
         prevBType = bonds[nActive - i];
         addEndAtom(thisPtr, prevPtr, prevBType, rosen_f, energy_f);
         #ifndef SIMP_NOPAIR
         system().pairPotential().addAtom(*thisPtr);
         #endif

         rosenbluth *= rosen_f;
         energy += energy_f;

         thisPtr += sign;
         prevPtr += sign;
      }

      // Step 2: grow the last atoms
      if (nActive >= 1) {
         prevPtr = thisPtr - sign;
         nextPtr = thisPtr + sign;
         prevBType = bonds[1];
         nextBType = bonds[0];

         // Invoke the orientation biased trimer re-bridging move
         addMiddleAtom(thisPtr, prevPtr, nextPtr,
                prevBType, nextBType, rosen_f, energy_f);
         #ifndef SIMP_NOPAIR
         system().pairPotential().addAtom(*thisPtr);
         #endif

         rosenbluth *= rosen_f;
         energy += energy_f;
      }

   }
   
   /* 
   * First determine the preferential orientation angle and angle biasing
   * spring constant; then compute the statistical weight for the biasing
   * function using Simpson's rule:
   *
   * Note: the variou branches represent certain geometrical extreme
   * situations.
   */
   void CfbRebridgeBase::orientationBiasTable
         (double l_20, double &prefAng, double &kappaAng, double &normConst)
   {
      double cos_ang;
   
      // Preferential orientation & effective spring constant
      if ( l_20 <= fabs(length21_ - length10_) ) {
   
         if (length21_ > length10_) {
            prefAng   = 0.0;
            kappaAng  = kappa10_ * (length21_ - l_20 - length10_);
            kappaAng *= length21_ * l_20 / (length21_ - l_20);
         } else if (length21_ == length10_) {
            prefAng   = acos(-1.0)/2.0;
            kappaAng  = 0.0;
         } else {
            prefAng   = acos(-1.0);
            kappaAng  = kappa10_ * (length10_ - l_20 - length21_);
            kappaAng *= length21_ * l_20 / (length21_ + l_20);
         }
   
      } else if ( l_20 >= (length21_ + length10_) ) { 
   
         prefAng  = 0.0;           
	 kappaAng = kappa10_ * (l_20 - length21_ - length10_) * length21_ * l_20
                     / (l_20 - length21_);	   
   
      } else {
   
         cos_ang  = (l_20*l_20 + length21_*length21_ - length10_*length10_)
                    * 0.5 / length21_ / l_20;
         prefAng  = acos(cos_ang);
         kappaAng = length21_ * l_20 / length10_;
         kappaAng = kappa10_ * kappaAng * kappaAng * (1.0 - cos_ang*cos_ang);
   
      }
   
      // Calculate normalization factor.
      if (kappaAng < 1.0E-7) {
         normConst = 2.0;
      } else {
         computeNormalization(prefAng, kappaAng, normConst);
      }
   }

   /* 
   * Use Simpson's rule to compute the normalization factor
   *   int { exp(-0.5*beta*kappaAng*(x-prefAng)**2)*sin(x) , {x, 0, pi} }
   */
   void CfbRebridgeBase::computeNormalization
      (double prefAng, double kappaAng, double &normConst)
   {
      int    i, N;                     // N must be even
      double delta, x, sum, tmp;
   
      // indices for Simpson's rule: 0,1,...,N
      N     = 400;                     // default value
      tmp   = 0.5*kappaAng;
      delta = acos(-1.0)/double(N);
   
      // the end points contribute nothing since sin(0)=sin(pi)=0.
      sum   = 0.0;
      x     = -delta;
      for (i=0; i < N-1; i=i+2) {
         x += delta;
         //sum += 2.0 * exp( -tmp * (x-prefAng)*(x-prefAng) ) * sin(x);
         sum += 2.0 * boltzmann( tmp * (x-prefAng)*(x-prefAng) ) * sin(x);
         x += delta;
         //sum += 4.0 * exp( -tmp * (x-prefAng)*(x-prefAng) ) * sin(x);
         sum += 4.0 * boltzmann( tmp * (x-prefAng)*(x-prefAng) ) * sin(x);
      }
      normConst = delta*sum/3.0;
   }

}
