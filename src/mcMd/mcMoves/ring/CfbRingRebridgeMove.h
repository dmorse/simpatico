#ifndef MCMD_CFB_RING_REBRIDGE_MOVE_H
#define MCMD_CFB_RING_REBRIDGE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/CfbRebridgeBase.h>  // base class
#include <util/containers/DArray.h>             // member template
#include <util/space/Vector.h>                   // member template parameter

namespace McMd
{

   using namespace Util;

   class McSystem;
   class Molecule;

   /**
   * Configuration bias rebridge moves for ring molecules.
   *
   * An interior bridge of atoms with length nRegrow is randomly selected and
   * is regrown using the configuration biased moves. The first nRegrow-1 atoms
   * are re-grown using the same biasing function as those for end atom moves,
   * and the last atom is re-grown using the crank-shaft like trimer re-growth
   * algorithm.
   *
   * \ingroup McMd_McMove_Module
   */
   class CfbRingRebridgeMove : public CfbRebridgeBase
   {
   
   public:
   
      /**
      * Constructor. 
      */
      CfbRingRebridgeMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

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

   private:

      /**
      * Ring version of deleting a sequence of atoms.
      * 
      * \param sign       direction to delete atoms.
      * \param molPtr     pointer to the active molecule.
      * \param endId      end atom ID of the active interior atom group.
      * \param bonds      bond types array.
      * \param rosenbluth rosenbluth factor.
      * \param energy     energy variation of the deleting operation.
      */
      void deleteSequence(int sign, Molecule *molPtr, int endId, int *bonds, 
                         double &rosenbluth, double &energy);

      /**
      * Ring version of growing a sequence of atoms.
      *
      * \param sign       Direction to delete atoms.
      * \param molPtr     Pointer to the active molecule.
      * \param beginId    Begin atom ID of the active interior atom group.
      * \param bonds      Bond types array.
      * \param rosenbluth Rosenbluth factor.
      * \param energy     Energy variation of the deleting operation.
      */
      void addSequence(int sign, Molecule *molPtr, int beginId, int *bonds, 
                         double &rosenbluth, double &energy);

      /**
      * Shift Atom index along a Ring.
      *
      * \param id Atom id to be regulated, which can be negative.
      * \param n  Number of atoms in the Ring.
      */
      int modId(int id, int n);

   };

   // Private inline methods.

   inline int CfbRingRebridgeMove::modId(int id, int n)
   {
      if (id >= n) id = id%n;
      if (id < 0) id = n - (-id)%n;
      return id;
   }

}      
#endif
