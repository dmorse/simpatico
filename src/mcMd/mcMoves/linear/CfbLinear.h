#ifndef MCMD_CFB_LINEAR_H
#define MCMD_CFB_LINEAR_H

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
   * Limitations: The current implementation works with model with 
   * external and angle potentials, but does not yet work with
   * dihedral potentials. The treatment of angle potentials does
   * not use the angle potential to bias the choice of trial bond
   * orientations, and so is not efficient for strong angle 
   * potentials. 
   *
   * \ingroup McMd_McMove_Module 
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
      * Read and validate speciesId and nTrial.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load and validate speciesId and nTrial.
      *
      * \param ar input (loading) archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save speciesId and nTrial.
      *
      * \param ar output (loading) archive
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * CFB algorithm for deleting an end atom from a Linear molecule.
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
      *   - energy is the total energy (bonded + nonbonded) of the end 
      *     atom before it was deleted.
      *
      * \param molecule  molecule
      * \param atomId  id of atom to be deleted
      * \param sign  end from which deletion is occuring (= +-1)
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
      * \param sign  end from which deletion is occuring (=+-1)
      * \param rosenbluth  Rosenbluth factor of added atom (out)
      * \param energy  potential energy of deleted atom (out)
      */
      void addAtom(Molecule& molecule, Atom& atom0, Atom& atom1, int atomId, 
                   int sign, double &rosenbluth, double &energy);

   protected:

      /**
      * Get nTrial.
      */
      int nTrial() const;

      /**
      * Get speciesId.
      */
      int speciesId() const;

   private:

      /// Maximum allowed number of trial positions
      static const int MaxTrial_ = 200;

      /// Integer identifier of molecular species
      int  speciesId_;

      /// Number of trial positions for regrown atom (<= MaxTrial_)
      int  nTrial_;

      #ifdef SIMP_ANGLE
      bool hasAngles_;
      #endif

      #ifdef SIMP_DIHEDRAL
      bool hasDihedrals_;
      #endif

      #ifdef SIMP_EXTERNAL
      bool hasExternal_;
      #endif

      /**
      * Validate and process parameters speciesId and nTrial.
      */
      void processParameters();

   };

   // Inline functions

   /*
   * Get nTrial.
   */
   inline int CfbLinear::nTrial() const
   {  return nTrial_; }

   /*
   * Get speciesId.
   */
   inline int CfbLinear::speciesId() const
   {  return speciesId_; }

}
#endif
