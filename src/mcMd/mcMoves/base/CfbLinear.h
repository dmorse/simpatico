#ifndef MCMD_CFB_END_BASE_H
#define MCMD_CFB_END_BASE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
   * This class provides methods deleteAtom and addAtom that
   * delete or add a single atom to a specified end of a linear chain.
   *
   * \ingroup McMd_McMove_Module MD_Module
   */
   class CfbLinear : public SystemMove
   {

   public:

      /**
      * Constructor.
      */
      CfbLinear(McSystem& system);

      /**
      * Destructor.
      */
      virtual ~CfbLinear();

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
      * \param molecule   molecule
      * \param atomId     id of atom to be deleted
      * \param sign       end from which deletion is occuring
      * \param rosenbluth nonbonded Rosenbluth factor of deleted atom (out)
      * \param energy     total potential energy of deleted atom (out)
      */
      void deleteAtom(Molecule& molecule, int atomId, int sign
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
      * \param molecule   molecule
      * \param atomId     id of atom to be deleted
      * \param sign       end from which deletion is occuring
      * \param rosenbluth Rosenbluth factor of added atom (out)
      * \param energy     potential energy of deleted atom (out)
      */
      void addAtom(Molecule& molecule, int atomId, int sign
                   double &rosenbluth, double &energy);

   protected:

      /// Maximum allowed number of trial positions for a regrown atom.
      static const int MaxTrial_ = 20;

      /// Actual number of trial positions for each regrown atom.
      int  nTrial_;

   };

}
#endif
