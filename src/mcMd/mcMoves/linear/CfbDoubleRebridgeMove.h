#ifndef MCMD_CFB_DOUBLE_REBRIDGE_MOVE_H
#define MCMD_CFB_DOUBLE_REBRIDGE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/CfbRebridgeBase.h> // base class
#include <util/containers/DArray.h>            // member template
#include <util/space/Vector.h>                  // member template parameter

namespace McMd
{

   using namespace Util;

   class McSystem;
   class Linear;

   /**
   * configuration bias trimer double rebridge moves, to reconnect two chains.
   *
   * \ingroup McMd_McMove_Module
   */
   class CfbDoubleRebridgeMove : public CfbRebridgeBase
   {
   
   public:
   
      /**
      * Constructor. 
      */
      CfbDoubleRebridgeMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
   
   protected:
   
      /// Array of old positions of temporarily deleted atoms (iMol).
      DArray<Vector> iOldPos_;

      /// Array of old positions of temporarily deleted atoms (jMol).
      DArray<Vector> jOldPos_;
    
      /// Integer index for molecular species.
      int speciesId_;
   
      /// Number of particles at end to attempt to regrow.
      int nRegrow_;

      /// upper bound for trial lenght of a bridge.
      double bridgeLength_;

   private:

      /// Scan potential bridging sites for old->new move
      bool forwardScan(int sign, int &iMol, int &jMol,
              int &beginId, double &prob);

      /// Scan potential bridging sites for new->old move
      bool reverseScan(int sign, int iMol, int jMol,
              int beginId, double &prob);
  
   };

}      
#endif
