#ifndef MCMD_RING_TETRA_REBRIDGE_MOVE_H
#define MCMD_RING_TETRA_REBRIDGE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/GroupRebridgeBase.h>  // base class
#include <util/space/Vector.h>                     // template parameter
#include <util/containers/DArray.h>               // member template

namespace McMd
{

   using namespace Util;

   class McSystem;
   class Molecule;

   /**
   * Molecule move that attempt the exchange of interior pieces of one close
   * rings. Analogous to the twist/untwist operations.
   *
   * \ingroup McMd_McMove_Module
   */
   class RingTetraRebridgeMove : public GroupRebridgeBase
   {
   
   public:
   
      /**
      * Constructor. 
      * 
      * \param system parent McSystem
      */
      RingTetraRebridgeMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
   
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive
      */
      virtual void save(Serializable::OArchive& ar);
  
   protected:
   
      /// Integer index for molecular species.
      int speciesId_;

      /// Number of atoms per molecule.
      int nAtom_;

      /// Upper bounds for trial length of a bridge.
      double upperBridge_;

      /// Lower bounds for trial length of a bridge.
      double lowerBridge_;

      /// Arrays of atom positions of the correct molecule shape.
      DArray<Vector> R_;

      /**
      * Scan potential rebridging sites.
      *
      * \param molPtr  Pointer of the molecule to be rebridged.
      * \param aId     The first atom Id of the first group to be rebridged.
      * \param cId     The first atom Id of the second group to be rebridged.
      */
      bool scanBridge(Molecule* molPtr, int aId, int &cId);

      /**
      * Shift Atom index along a Ring.
      *
      * \param id Atom id to be regulated, which can be negative.
      * \param n  Number of atoms in the Ring.
      */
      int modId(int id, int n);

   };

   // Private inline methods.

   inline int RingTetraRebridgeMove::modId(int id, int n)
   {
      if (id >= n) id = id%n;
      if (id < 0) id = n - (-id)%n;
      return id;
   }

}      
#endif
