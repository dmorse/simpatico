#ifndef DDMD_OUTPUT_ENERGY_H
#define DDMD_OUTPUT_ENERGY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically write simulation energies to file.
   *
   * \sa \ref ddMd_analyzer_OutputEnergy_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class OutputEnergy : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      OutputEnergy(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~OutputEnergy()
      {} 
   
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
      * Clear nSample counter.
      */
      virtual void clear();
  
      /**
      * Setup - open output file.
      */
      virtual void setup();

      /**
      * Write energy to file
      *
      * \param iStep MC step index
      */
      virtual void sample(long iStep);

   private:
 
      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long    nSample_;
   
      /// Has readParam been called?
      long    isInitialized_;
   
   };

}
#endif 
