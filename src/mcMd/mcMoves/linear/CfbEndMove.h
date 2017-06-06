#ifndef MCMD_CFB_END_MOVE_H
#define MCMD_CFB_END_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/CfbEndBase.h>  // base class
#include <util/containers/DArray.h>        // member template
#include <util/space/Vector.h>              // member template parameter

namespace Simp {
   class Linear;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   class McSystem;

   /**
   * Configuration bias end regrowth move for flexible linear chains.
   *
   * \ingroup McMd_McMove_Module
   */
   class CfbEndMove : public CfbEndBase
   {
   
   public:
   
      /**
      * Constructor. 
      */
      CfbEndMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
   
   protected:
   
      /// Array of old positions of temporarily deleted atoms (temporary).
      DArray<Vector> oldPos_;
    
      /// Integer index for molecular species.
      int speciesId_;
   
      /// Number of particles at end to attempt to regrow
      int nRegrow_;
   
      /// Make unit test class a friend
      friend class CbEndMoveTest;
   
   };

}      
#endif
