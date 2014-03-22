#ifndef MCMD_VIRIAL_STRESSTENSOR_AVERAGE_H
#define MCMD_VIRIAL_STRESSTENSOR_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/Analyzer.h>
#include <mcMd/analyzers/SystemAnalyzer.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <util/space/Tensor.h>
#include <util/accumulators/Average.h>

namespace McMd
{

   using namespace Util;

   /**
   * Periodically write (tensor) StressTensor to file.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class McVirialStressTensorAverage : public SystemAnalyzer<McSystem>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent McSystem object. 
      */
      McVirialStressTensorAverage(McSystem& system);
   
      /**
      * Destructor.
      */
      virtual ~McVirialStressTensorAverage()
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
      * Sample virial stress to accumulators
      *
      * \param iStep MD step index
      */
      virtual void sample(long iStep);

      /**
      * Dump configuration to file
      */
      virtual void output();

   private:
 
      /// Output file stream
      std::ofstream outputFile_;
      
      /// Average object to save sxx
      Average sxxAccumulator_;

      /// Average object to save sxx
      Average sxyAccumulator_;

      /// Average object to save sxx
      Average sxzAccumulator_;

      /// Average object to save sxx
      Average syxAccumulator_;

      /// Average object to save sxx
      Average syyAccumulator_;

      /// Average object to save sxx
      Average syzAccumulator_;

      /// Average object to save sxx
      Average szxAccumulator_;

      /// Average object to save sxx
      Average szyAccumulator_;

      /// Average object to save sxx
      Average szzAccumulator_;

      /// Number of samples per block average output
      int nSamplePerBlock_;

      /// Number of configurations dumped thus far (first dump is zero).
      long    nSample_;
   
      /// Has readParam been called?
      long    isInitialized_;
   
   };

}
#endif 
