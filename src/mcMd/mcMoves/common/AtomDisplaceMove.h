#ifndef MCMD_ATOM_DISPLACE_MOVE_H
#define MCMD_ATOM_DISPLACE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>        // base class
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * Random displacement of one atom.
   *
   * \sa \ref mcMd_mcMove_AtomDisplaceMove_page "parameter file format"
   *
   * \ingroup McMd_McMove_Module McMove_Module
   */
   class AtomDisplaceMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      */
      AtomDisplaceMove(McSystem& system);
   
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
      * Generate, attempt and accept or reject a move.
      */
      virtual bool move();

   private:

      /// Maximum magnitude of displacement.
      double delta_;

      /// Integer Id of Species.
      int    speciesId_;
   
   };

}      
#endif
