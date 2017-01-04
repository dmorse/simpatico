#ifndef MCMD_RIGID_DISPLACE_MOVE_H
#define MCMD_RIGID_DISPLACE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>  // base class
#include <util/containers/DArray.h>   // member
#include <util/space/Vector.h>         // member template parameter
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * Random rigid displacement of a molecule.
   *
   * \sa \ref mcMd_mcMove_RigidDisplaceMove_page "parameter file format"
   *
   * \ingroup McMd_McMove_Module McMove_Module
   */
   class RigidDisplaceMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      */
      RigidDisplaceMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      *
      * \param in input parameter stream
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
      * Serialize internal state to/from archive.
      *
      * \param ar       archive
      * \param version  id for archive version
      */
      template <class Archive>
      void serialize(Archive &ar, const unsigned int version);

      /**
      * Generate, attempt and accept or reject a move.
      */
      virtual bool move();

   private:

      /// Array of old positions.
      DArray<Vector> oldPositions_;

      /// Maximum magnitude of displacement
      double delta_;

      /// Integer Id of Species.
      int    speciesId_;
   
      /// Number of atoms per molecule for this species.
      int    nAtom_;
   
   };

}      
#endif
