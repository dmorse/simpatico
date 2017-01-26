#ifndef MCMD_CFB_END_BASE_H
#define MCMD_CFB_END_BASE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>  // base class

namespace McMd
{

   using namespace Util;

   class Atom;
   class McSystem;

   /**
   * Base class for configuration bias (CFB) end regrowth moves.
   *
   * This class provides methods deleteEndAtom and addEndAtom that
   * delete or add a single atom to a specified end of a linear chain.
   * These are the basic building blocks of CFB moves for flexible
   * linear chains. 
   *
   * \ingroup McMd_McMove_Module MD_Module
   */
   class CfbEndBase : public SystemMove 
   {

   public:

      /**
      * Constructor. 
      */
      CfbEndBase(McSystem& system);
   
      /**
      * Destructor. 
      */
      virtual ~CfbEndBase();

      /**
      * Read parameter nTrial.
      *
      * This function is used only for testing. Subclasses read the
      * protected member nTrial directly. 
      */
      virtual void readParameters(std::istream& in);

      // No loadParameters or save methods are needed. The nTrial
      // can be directly loaded and saved by subclasses.
   
      /**
      * CFB algorithm for deleting an end atom from a flexible chain.
      *
      * This function computes the energy of an end atom and a Rosenbluth 
      * factor for removing it. It does not remove the end atom from the 
      * system cell list.
      *
      * Upon return:
      *   
      *   - rosenbluth is the nonbonded Rosenblush factor for the deleted
      *     atom, i.e., the sum of Boltzmann factors from nonbonded pair
      *     interactions for the initial position and nTrial_ - 1 trials.
      *
      *   - energy is the total energy (bonded + nonbonded) of the end atom 
      *     before it was deleted.
      *
      * \param endPtr     ptr to end atom, which we attempt to remove
      * \param pvtPtr     ptr to atom next to end, or end after removal
      * \param bondType   type of bond connecting pvt and end atoms
      * \param rosenbluth nonbonded Rosenbluth factor of deleted atom (out)
      * \param energy     total potential energy of deleted atom (out)
      */
      void deleteEndAtom(Atom* endPtr, Atom* pvtPtr, int bondType,
                         double &rosenbluth, double &energy);
   
  
      /**
      * Configuration bias algorithm for adding an atom to a chain end.
      *
      * This function generates and computes Rosenbluth factors for nTrial 
      * trial positions, chooses one, updates the atomic position. It does
      * not add the end atom to the system cell list.
      *
      * Upon return:
      *  
      *   - rosenbluth is the nonbonded Rosenblush factor for the added 
      *     atom, i.e., the sum of Boltzmann factors from nonbonded pair 
      *     interactions for all nTrial_ trial positions.
      *
      *   - energy is the total energy (bonded + nonbonded) of the new end
      *     atom in its chosen position.
      *
      * \param endPtr     ptr to new end atom, which we attempt to add
      * \param pvtPtr     end atom of current chain, next to end
      * \param bondType   type of bond connecting pvt and end atoms
      * \param rosenbluth Rosenbluth factor of added atom (out)
      * \param energy     potential energy of deleted atom (out)
      */
      void addEndAtom(Atom* endPtr, Atom* pvtPtr, int bondType,
                      double &rosenbluth, double &energy);
   
   protected:

      /// Maximum allowed number of trial positions for a regrown atom.
      static const int MaxTrial_ = 20; 
   
      /// Actual number of trial positions for each regrown atom.
      int  nTrial_; 
   
      /// Grant friend access to unit test class
      //  friend class CbEndBaseTest;

   };

}      
#endif
