#ifndef  SIMP_NOPAIR
#ifndef MCMD_MD_PAIR_ENERGY_COEFFICIENTS_H
#define MCMD_MD_PAIR_ENERGY_COEFFICIENTS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/mdSimulation/MdSystem.h>         // base template parameter
#include <mcMd/analyzers/util/PairSelector.h>

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
   * Analyzer to output the total pair energy and the sum
   * of squares of the monomeric and molecular pair energy
   *
   * See \ref mcMd_analyzer_MdPairEnergyCoefficients_page "here" 
   * for the parameter file format and any other user documentation.
   * 
   * \ingroup McMd_Analyzer_Md_Module
   */
   class MdPairEnergyCoefficients : public SystemAnalyzer<MdSystem>
   {

   public:
   
      /// Constructor.
      MdPairEnergyCoefficients(MdSystem& system);

      /// Destructor
      ~MdPairEnergyCoefficients();

      /**
      * Read parameters and initialize.
      *
      * Reads output file, pair selector and maximum number of neighbors
      * per molecule
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);
 
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /// Evaluate energy and print.
      virtual void sample(long iStep);

      /// Output final summary and file format
      virtual void output();

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

      /// Array of neighbor lists for every molecule in every species
      DArray< DArray< DSArray<  Pair< Atom *> > > > moleculeNeighbors_;

      /// Pair energy with every neighboring molecule in every species
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
      ar & pairEnergyAccumulator_;
      ar & moleculePESqAccumulator_;
      ar & twoMoleculePESqAccumulator_;
      ar & pESqAccumulator_;
      //DArray< DArray< DSArray<  Pair< Atom *> > > > moleculeNeighbors_;
      //DArray< DArray< double > > twoMoleculePairEnergy_;

   }

}
#endif 
#endif
