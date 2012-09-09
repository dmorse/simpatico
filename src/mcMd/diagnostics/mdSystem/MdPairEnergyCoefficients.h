#ifndef  INTER_NOPAIR
#ifndef MCMD_MD_PAIR_ENERGY_COEFFICIENTS_H
#define MCMD_MD_PAIR_ENERGY_COEFFICIENTS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/mdSimulation/MdSystem.h>         // base template parameter
#include <mcMd/diagnostics/util/PairSelector.h>

#include <util/global.h> 
#include <util/containers/Pair.h>
#include <util/containers/DSArray.h>
#include <util/accumulators/Average.h>

namespace McMd
{

   using namespace Util;

   class MdPairPotential;
   class PairList;

   /**
   * Diagnostic to output the total pair energy and the
   * sum of squares of the monomeric and molecular pair energy
   * 
   * \ingroup McMd_Diagnostic_Module
   */
   class MdPairEnergyCoefficients : public SystemDiagnostic<MdSystem>
   {

   public:
   
      /// Constructor.
      MdPairEnergyCoefficients(MdSystem& system);

      /// Destructor
      ~MdPairEnergyCoefficients();

      /// Read output file, pair selector and maximum number of neighbors
      /// per molecule
      virtual void readParameters(std::istream& in);
 
      /// Evaluate energy and print.
      virtual void sample(long iStep);

      /// Output final summary and file format
      virtual void output();

      /**
      * Save state to binary file archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchiveType& ar);

      /**
      * Load state from a binary file archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void load(Serializable::IArchiveType& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// Rule specifying which atom pairs to accept
      PairSelector  selector_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int nAtomType_;

      /// Number of species, copied from Simulation::nSpecices().
      int nSpecies_;

      /// Pointer to the pair list of the associated MdSystem
      const PairList  *pairListPtr_;

      /// Pointer to the pair potential of the associated System
      MdPairPotential* pairPotentialPtr_;

      /// Pointer to the boundary of the system
      Boundary *boundaryPtr_;

      /// Maximum number of neighbors per molecule
      int maxMoleculeNeighbors_;

      /// array of neighbor lists for every molecule in every species
      DArray< DArray< DSArray<  Pair< Atom *> > > > moleculeNeighbors_;

      /// pair energy with every neighboring molecule in every species
      /// for a given molecule
      DArray< DArray< double > > twoMoleculePairEnergy_;

      /// Accumulator for monomer pair energy
      Average pairEnergyAccumulator_;

      /// Accumulator for squares of molecular pair energy
      Average moleculePESqAccumulator_;

      /// Accumulator for squares of two molecule pair energy
      Average twoMoleculePESqAccumulator_;

      /// Accumulator for variance of total pair energy
      Average pESqAccumulator_;

      /// Output file stream
      std::ofstream outputFile_;

      // Has readParam been called?
      bool isInitialized_;

      /*
      * Methods
      */

      // Reset the internal molecules' neighbor list arrays
      void clear();

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void MdPairEnergyCoefficients::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }

      ar & pairEnergyAccumulator_;
      ar & moleculePESqAccumulator_;
      ar & twoMoleculePESqAccumulator_;
      ar & pESqAccumulator_;
      serializeCheck(ar, nAtomType_, "nAtomType");
      serializeCheck(ar, nSpecies_, "nSpecies");
      serializeCheck(ar, maxMoleculeNeighbors_, "maxMoleculeNeighbors");
      //PairSelector  selector_;
      //DArray< DArray< DSArray<  Pair< Atom *> > > > moleculeNeighbors_;
      //DArray< DArray< double > > twoMoleculePairEnergy_;

   }

}
#endif 
#endif
