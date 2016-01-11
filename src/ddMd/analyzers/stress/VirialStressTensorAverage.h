#ifndef DDMD_VIRIAL_STRESS_TENSOR_AVERAGE_H
#define DDMD_VIRIAL_STRESS_TENSOR_AVERAGE_H

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
#include <util/accumulators/Average.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Compute average and output block averages of virial stress tensor.
   *
   * \sa \ref ddMd_analyzer_VirialStressTensorAverage_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Stress_Module
   */
   class VirialStressTensorAverage : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      VirialStressTensorAverage(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~VirialStressTensorAverage()
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
