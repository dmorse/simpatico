#ifndef MCMD_CFB_REBRIDGE_MOVE_H
#define MCMD_CFB_REBRIDGE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/CfbRebridgeBase.h>  // base class
#include <util/containers/DArray.h>             // member template
#include <util/space/Vector.h>                  // member template parameter

namespace Simp {
   class Linear;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   class McSystem;

   /**
   * Config-bias move for internal segment of a flexible linear polymer.
   *
   * \sa \ref mcMd_mcMove_CfbRebridgeMove_page "parameter file format"
   *
   * \ingroup McMd_McMove_Module
   */
   class CfbRebridgeMove : public CfbRebridgeBase
   {
   
   public:
   
      /**
      * Constructor. 
      */
      CfbRebridgeMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
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
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
   
   protected:
   
      /// Array of old positions of temporarily deleted atoms
      DArray<Vector> oldPos_;
    
      /// Integer index for molecular species.
      int speciesId_;
   
      /// Number of particles at end to attempt to regrow
      int nRegrow_;
   
   };

}      
#endif
