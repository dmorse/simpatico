#ifndef SIMP_SPECIES_FINDER_H
#define SIMP_SPECIES_FINDER_H
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2020, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>

namespace Simp
{

   using namespace Util;

   /**
   * Find the context of molecule, atom or group from a global object id.
   *
   * This class is designed to compute the molecular context of a molecule
   * or a part of a molecule (i.e., an atom or covalent group), given the
   * value of a global id for the molecule or part that is known to have
   * been generated as part of an ordered sequence in which molecules 
   * of the same species are listed sequentially, and parts within each
   * molecule listed sequentially. The class must be initialized by
   * providing the number of species, the number of molecules in each
   * species and the number of parts (atoms or groups) within each 
   * molecule of a each species. 
   *
   * Initialization:
   * \code
   *   
   *   SpeciesFinder finder; 
   *   finder.allocate(nSpecies);
   *   for (int i=0; i < nSpecies; ++i) {
   *      finder.setSpecies(nMolecule, nPart);
   *   }
   *   finder.initialize();
   * 
   * \endcode
   *
   * The functions findMolecule and findPart may only be called after
   * this initialization process.
   *
   * \ingroup Simp_Species_Module
   */
   class SpeciesFinder 
   {

   public:

      /**
      * Context (local indices) for a molecule.
      */
      struct Molecule {

         /// Index of parent molecule species
         int speciesId;

         /// Index of molecule within its species, from 0.
         int moleculeId;

      };

      /**
      * Context (local indices) for a molecule part (atom or group).
      */
      struct Part {

         /// Index of parent molecule species
         int speciesId;

         /// Index of molecule within its species, from 0.
         int moleculeId;

         /// Index of atom or group within its molecule, from 0.
         int partId;

      };

      /**
      * Constructor.
      */
      SpeciesFinder();

      /**
      * Allocate all required private arrays.
      *
      * \param nSpecies number of molecular species
      */
      void allocate(int nSpecies);

      /**
      * Set data for one species.
      *
      * \param iSpecies species index
      * \param nMolecule number of molecules of this species
      * \param nPart number of parts (atoms or groups) per molecule
      */
      void setSpecies(int iSpecies, int nMolecule, int nPart);

      /**
      * Finish initialization procedure, after adding all species.
      */
      void initialize();

      /**
      * Find the context of a molecule.
      *
      * \param moleculeId Global id for a molecule (input)
      * \param context molecule context (output)
      */
      void findMolecule(int moleculeId, Molecule& context);

      /**
      * Find the context of a part of a molecule (atom or group)
      *
      * \param moleculeId Global id for the part (input)
      * \param context atom or group context (output)
      */
      void findPart(int partId, Part& context);

   private:

      /// Number of species
      int nSpecies_;

      /// Number of molecules of each species.
      DArray<int> nMolecule_;

      // Number of parts (atoms or groups) in each molecule.
      DArray<int> nPart_;

      /// First molecule id of each species.
      DArray<int> firstMoleculeId_;

      /// First part id of first molecule of each species. 
      DArray<int> firstPartId_;

   };

} 
#endif
