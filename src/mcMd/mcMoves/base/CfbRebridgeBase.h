#ifndef MCMD_CFB_REBRIDGE_BASE_H
#define MCMD_CFB_REBRIDGE_BASE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbEndBase.h"        // base class
#include <util/space/Vector.h>  // math utility class

namespace McMd
{

   using namespace Util;

   class McSystem;
   class Atom;

   /**
   * Base class configuration bias moves internal segment regrowth moves.
   *
   * Base class for configuration bias Monte Carlo moves that erase and
   * regrow internal segments of linear polymer chains, with fixed
   * endpoints.
   *
   * \ingroup McMd_McMove_Module
   */
   class CfbRebridgeBase : public CfbEndBase 
   {
   
   public:
   
      /**
      * Constructor. 
      */
      CfbRebridgeBase(McSystem& system);
   
      /**
      * Destructor. 
      */
      virtual ~CfbRebridgeBase();

      /**
      * Read species type, nTrial, and parameters needed to evaluate the
      * orientation bias.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);
   
      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Initialize the arrays for preferential angle, effective spring constant
      * and normalization factor.
      */
      virtual void setup();
   
   protected:

      /// Equilibrium bond lengths (2-1 bond).
      double length21_;

      /// Equilibrium bond lengths (1-0 bonds). 
      double length10_;

      /// Spring constant (1-0 bond).
      double kappa10_;

      /**
      * Configuration bias algorithm for deleting a particle between two
      * monomers and compute the corresponding Rosenbluth factor.
      */
      void deleteMiddleAtom(Atom* partPtr, Atom* prevPtr, Atom* nextPtr,
             int prevBType, int nextBType, double &rosenbluth, double &energy);

      /**
      * Configuration bias algorithm for adding the last particle with special
      * bias in favor of closing the bridge.
      */
      void addMiddleAtom(Atom* partPtr, Atom* prevPtr, Atom* nextPtr,
             int prevBType, int nextBType, double &rosenbluth, double &energy);

      /**
      * Configuration bias algorithm for deleting a consequtive sequence of atoms.
      */
      void deleteSequence(int nActive, int sign, Atom* endPtr, int *bonds,
                          double &rosenbluth, double &energy);
   
      /**
      * Configuration bias algorithm for adding a consequtive sequence of atoms.
      */
      void addSequence(int nActive, int sign, Atom* beginPtr, int *bonds,
                          double &rosenbluth, double &energy);

   private:
   
      /**
      * Pre-computed preferential angle, orientation bias spring constant and
      * the normalization factor for the angular biasing function. The values
      * are tabulated for various values of l_20.
      */

      /// Maximum number of l_20 values in the table (>3).
      static const int MaxBin_ = 102;

      /// A tabulated l20 values.
      double l20Table[MaxBin_];

      /// Preferential angle.
      double angTable[MaxBin_];

      /// Effective spring constant.
      double kappaTable[MaxBin_];

      /// Normalization factor.
      double normalTable[MaxBin_];
   
      /**
      * Inline method: get the preferential angle and effective spring constant
      * for given l_20.
      */
      inline void getAngKappaNorm(double l_20, double &preAng,
                                  double &kappa, double &normConst);
  
      /**
      * Inline method: Compute the orientation bias factor for Crank-Shaft type
      * move based on orientation angle.
      */
      inline void orientationBias(Vector u_21, Vector u_20,
          double prefAng, double kappaAng, double &bias);   

      /**
      * Compute preferential orientation angle, angle biasing spring constant
      * and the normalization constant, based on the bonding length between
      * particle 2 and 0.
      */
      void orientationBiasTable
         (double l_20, double &prefAng, double &kappaAng, double &normConst);
   
      /**
      * Compute the normalization factor for given preferential angle and
      * effective spring constant.
      */
      void computeNormalization
         (double prefAng, double kappaAng, double &normConst);
 
   };

   /* 
   * Inline method: get the preferential angle and effective spring constant for
   * given l_20
   */
   void CfbRebridgeBase::getAngKappaNorm(double l_20, double &prefAng,
           double &kappa, double &normConst)
   {
      if (l_20 <= l20Table[0]) {
         prefAng   = angTable[0];
         kappa     = kappaTable[0];
         normConst = normalTable[0];
      } else if (l_20 >= l20Table[MaxBin_-1]) {
         prefAng   = angTable[MaxBin_-1];
         kappa     = kappaTable[MaxBin_-1];
         normConst = normalTable[MaxBin_-1];
      } else {
         int i;
         i = int ( (l_20 - l20Table[0]) / (l20Table[2] - l20Table[1]) );
         prefAng   = angTable[i+1];
         kappa     = kappaTable[i+1];
         normConst = normalTable[i+1];
      }
   }
   
   /* 
   * Inline method: compute the actual angle and the Boltzmann weight.
   */
   void CfbRebridgeBase::orientationBias(
      Vector u_21, Vector u_20, double prefAng, double kappaAng, double &bias)
   {
      double actualAng;
   
      actualAng = acos(u_21.dot(u_20));
      bias = boltzmann(0.5*kappaAng*(actualAng - prefAng)*(actualAng - prefAng));
   }
   
}
#endif
