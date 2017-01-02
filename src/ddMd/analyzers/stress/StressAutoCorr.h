#ifndef DDMD_STRESS_AUTO_CORR_H
#define DDMD_STRESS_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/simulation/Simulation.h>
#include <util/mpi/MpiLoader.h>
#include <util/space/Tensor.h>
#include <util/accumulators/AutoCorrelation.h>     // member template

namespace DdMd
{

   using namespace Util;

   /**
   * Compute stress autocorrelation function for a liquid.
   *
   * \ingroup DdMd_Analyzer_Stress_Module
   */
   class StressAutoCorr : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      StressAutoCorr(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~StressAutoCorr()
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
      * Setup accumulator!
      */
      virtual void setup();
  
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
      std::ofstream  outputFile_;
      
      /// Statistical accumulator.
      AutoCorrelation<Tensor, double>*  accumulatorPtr_;

      /// Buffer capacity per stage (# values stored)
      int  bufferCapacity_;

      /// Maximum stage index for descendant AutoCorrStage objects
      int  maxStageId_;

      /// Has readParam been called?
      long  isInitialized_;
   
   };

}
#endif 
