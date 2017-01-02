#ifndef DDMD_PAIR_ENERGY_ANALYZER_H
#define DDMD_PAIR_ENERGY_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/AverageAnalyzer.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Compute and analyze pair energies.
   *
   * \sa \ref ddMd_analyzer_PairEnergyAnalyzer_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Energy_Module
   */
   class PairEnergyAnalyzer : public AverageAnalyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      PairEnergyAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~PairEnergyAnalyzer(); 
   
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
  
   protected:

      /**
      * Function to compute value.
      *
      * Call on all processors.
      */
      virtual void compute();

      /**
      * Current value, set by compute function.
      *
      * Call only on master.
      */
      virtual double value();

   private:
   
      /** 
      * Pair of atom type ids.
      */
      FArray<int, 2>  typeIdPair_;

   };

}
#endif 
