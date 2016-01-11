#ifndef DDMD_PAIR_ENERGY_SELECTOR_H
#define DDMD_PAIR_ENERGY_SELECTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/analyzers/util/PairSelector.h>
#include <ddMd/neighbor/PairList.h>
#include <ddMd/neighbor/CellList.h>
#include <ddMd/neighbor/PairIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/simulation/Simulation.h>
#include <util/accumulators/Average.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically write simulation energies to Log output.
   *
   * \sa \ref ddMd_analyzer_PairEnergySelector_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class PairEnergySelector : public Analyzer
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      PairEnergySelector(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~PairEnergySelector();

      /**
      * Read dumpPrefix and interval.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Dump configuration to file
      *
      * \param iStep MD step index
      */
      virtual void sample(long iStep);

      /**
      * Dump configuration to file
      */
      virtual void output();

   private:

      /// Output file stream.
      std::ofstream  outputIntra_;

      /// Output file stream.
      std::ofstream  outputInter_;

      /// Pointer to Pair Potential
      PairPotential  *pairPotentialPtr_;

      /// Pointer to Pair List
      PairList       *pairListPtr_;

      /// Number of atom types
      int nAtomType_;

      /// Pair selectors
      DMatrix<PairSelector>   selectors_;

      /// Array of pairEnergies for each selector
      DArray<Average>   pairEnergies_;

      /// Has readParam been called?
      bool  isInitialized_;

   };

}
#endif 
