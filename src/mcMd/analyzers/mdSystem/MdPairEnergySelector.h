#ifndef  INTER_NOPAIR
#ifndef MCMD_MD_PAIR_ENERGY_SELECTOR_H
#define MCMD_MD_PAIR_ENERGY_SELECTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
   * Analyzer to output the total intermolecular pair energy
   * for the purpose of calculating the intermolecular correlation
   * number
   * 
   * \ingroup McMd_Analyzer_Md_Module
   */
   class MdPairEnergySelector : public SystemAnalyzer<MdSystem>
   {

   public:

      /// Constructor.
      MdPairEnergySelector(MdSystem& system);

      /// Destructor
      ~MdPairEnergySelector();

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

      /// Number of pair types
      int nPairTypes_;

      /// Array of PairSelector objects to calculate NerG n stuff
      DArray<PairSelector> selectors_;

      /// Array of pairEnergies for each selector
      DArray<Average>   pairEnergies_;

      /// Pointer to the pair list of the associated MdSystem
      const PairList  *pairListPtr_;

      /// Pointer to the pair potential of the associated System
      MdPairPotential* pairPotentialPtr_;

      /// Pointer to the boundary of the system
      Boundary *boundaryPtr_;

      /// Output file stream
      std::ofstream outputFile_;

      // Has readParam been called?
      bool isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void MdPairEnergySelector::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nPairTypes_;
//      ar & selectors_;

      int i;
      for (i = 0; i < nPairTypes_; i++) {
         ar & selectors_[i];
      }
      for (i = 0; i <= nPairTypes_; i++) {
         ar & pairEnergies_[i];
      }

//      ar & pairEnergies_;
      //DArray< DArray< DSArray<  Pair< Atom *> > > > moleculeNeighbors_;
      //DArray< DArray< double > > twoMoleculePairEnergy_;

   }

}
#endif 
#endif
