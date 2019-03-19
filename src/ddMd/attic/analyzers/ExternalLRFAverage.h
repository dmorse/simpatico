#ifndef DDMD_EXTERNAL_LRF_AVERAGE_H
#define DDMD_EXTERNAL_LRF_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/simulation/Simulation.h>
#include <util/accumulators/Average.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically write simulation energies to Log output.
   *
   * \sa \ref ddMd_analyzer_ExternalLRFAverage_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class ExternalLRFAverage : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      ExternalLRFAverage(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~ExternalLRFAverage(); 
   
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
      std::ofstream  outputFile_;

      /// Average object is to be set in the master processor!
      Average  *accumulator_;

      /// Number of samples per block average output.
      int  nSamplePerBlock_;
   
      /// Has readParam been called?
      bool  isInitialized_;
   
   };

}
#endif 
