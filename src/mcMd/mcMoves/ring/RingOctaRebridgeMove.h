#ifndef MCMD_RING_OCTA_REBRIDGE_MOVE_H
#define MCMD_RING_OCTA_REBRIDGE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/GroupRebridgeBase.h>  // base class
#include <mcMd/neighbor/CellList.h>               // member variable

namespace McMd
{

   using namespace Util;

   class McSystem;
   class Molecule;

   /**
   * Molecule rebridging move for a Ring species.
   *
   * The exchange of positions of two closely approaching atoms is attempted.
   * Each atom is bonded to two other atoms. The move requires that the 6 atoms
   * form an octahedron in which the edge lengths satisfying certain distance
   * criterions, before the exchange is attempted.
   *
   * \ingroup McMd_McMove_Module
   */
   class RingOctaRebridgeMove : public GroupRebridgeBase
   {
   
   public:
   
      /**
      * Constructor. 
      */
      RingOctaRebridgeMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
   
   protected:
   
      /// Integer index for molecular species.
      int speciesId_;
   
      /// Upper bounds for trial length of a bridge.
      double upperBridge_;

      /// Lower bounds for trial length of a bridge.
      double lowerBridge_;

      /// Neighbor list around test atom.
      CellList::NeighborArray neighbors_;

      /**
      * Scan potential rebridging sites
      *
      * \param molPtr  Pointer to the first molecule.
      * \param mId     The middle atom Id on the first molecule.
      * \param molId2  The index of the second molecule.
      * \param nId     The middle atom Id on the second molecule.
      */
      bool scanBridge(Molecule* molPtr, int mId, int &molId2, int &nId);

      /**
      * Shift Atom index along a Ring.
      *
      * \param id  Atom id to be regulated, which can be negative.
      * \param n   Number of atoms in the Ring.
      */
      int modId(int id, int n);

   };

   // Private inline methods.

   inline int RingOctaRebridgeMove::modId(int id, int n)
   {
      if (id >= n) id = id%n;
      if (id < 0) id = n - (-id)%n;
      return id;
   }

}      
#endif
