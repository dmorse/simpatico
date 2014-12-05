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
   class Molecule;
   class McSystem;

   /**
   * Base class for configuration bias (CFB) end regrowth moves.
   *
   * This class provides functions that use Rosenbluth sampling to 
   * delete or add a single atom to a specified end of a Linear 
   * molecule. These are the building blocks of configuration bias
   * regrowth and reptation algorithms for linear chains.
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
      * CFB algorithm for deleting an end atom from a Linear molecule.
      *
      * This function computes the energy of an end atom and a Rosenbluth
      * factor for removing it. It does not remove the end atom from the
      * system cell list.
      *
      * Note: Reference to atom0 and atom1 are passed as parameter to allow 
      * the use of this algorithm in reptation algorithm, in which a memory 
      * location from the tail (removal) end is used to temporarily 
      * represented the head (added) atom.
      *
      * Upon return:
      *
      *   - rosenbluth is the nonbonded Rosenblush factor for the deleted
      *     atom, i.e., the sum of Boltzmann factors from nonbonded pair
      *     interactions for the initial position and nTrial_ - 1 trials.
      *
      *   - energy is the total energy (bonded + nonbonded) of the end 
      *     atom before it was deleted.
      *
      * \param molecule  molecule
      * \param atomId  id of atom to be deleted
      * \param sign  end from which deletion is occuring
      * \param rosenbluth  nonbonded Rosenbluth factor of deleted atom (out)
      * \param energy  total potential energy of deleted atom (out)
      */
      void deleteAtom(Molecule& molecule, int atomId, int sign, 
                      double &rosenbluth, double &energy);

      /**
      * Configuration bias algorithm for adding an atom to a Linear molecule.
      *
      * This function generates and computes Rosenbluth factors for nTrial
      * trial positions, chooses one, updates the atomic position. It does
      * not add the end atom to the system cell list.
      *
      * The function signature requires references to atom0 (the atom to be
      * added) and atom1 (the one to which it is bonded) in order to allow
      * this to be used in a reptation algorithm in which atoms are not 
      * added in place: In the reptation algorithm, a trial removal is made
      * at the tail end of the molecule, and then the removed atom is used
      * to store the position of trial addition at the head end. The atom
      * positions are then shifted to restore the proper order in memory
      * only if the move is accepted. 
      *
      * Upon return:
      *
      *   - rosenbluth is assigned the nonbonded Rosenblush factor for the 
      *     added atom, i.e., the sum of Boltzmann factors for all nTrial_
      *     trial positions.
      *
      *   - energy is assigned the total energy (bonded + nonbonded) of 
      *     the new atom in the selected new position.
      *
      * \param molecule  molecule
      * \param atom0  atom to be added
      * \param atom1  atom to which atom0 is bonded (pivot atom)
      * \param atomId  local id of atom0 (0 <= id < molecule.nAtom())
      * \param sign  end from which deletion is occuring
      * \param rosenbluth  Rosenbluth factor of added atom (out)
      * \param energy  potential energy of deleted atom (out)
      */
      void addAtom(Molecule& molecule, Atom& atom0, Atom& atom1, int atomId, 
                   int sign, double &rosenbluth, double &energy);

   protected:

      /**
      * Read and validate speciesId parameter.
      *
      * \param in input parameter file
      */
      void readSpeciesId(std::istream& in);

      /**
      * Read and validate nTrial parameter.
      *
      * \param in input parameter file
      */
      void readnTrial(std::istream& in);

      /**
      * Get nTrial.
      */
      int nTrial() const
      {  return nTrial_; }

      /**
      * Get speciesId.
      */
      int nspeciesId() const
      {  return nTrial_; }

   private:

      /// Maximum allowed number of trial positions for a regrown atom.
      static const int MaxTrial_ = 20;

      /// Integer identifier of molecular species
      int  speciesId_;

      /// Actual number of trial positions for each regrown atom.
      int  nTrial_;

      #ifdef INTER_ANGLE
      bool hasAngles_;
      #endif

      #ifdef INTER_DIHEDRAL
      bool hasDihedrals_;
      #endif

      #ifdef INTER_EXTERNAL
      bool hasExternal_;
      #endif

   };

}
#endif
