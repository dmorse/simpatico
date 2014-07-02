#ifndef DDMD_VIRIAL_STRESS_TENSOR_H
#define DDMD_VIRIAL_STRESS_TENSOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * Periodically write (tensor) StressTensor to file.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class VirialStressTensor : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      VirialStressTensor(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~VirialStressTensor()
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

   private:

      ///Output file stream
      std::ofstream outputFile_;
   
      /// Number of samples
      int   nSample_;

      /// Has readParam been called?
      bool  isInitialized_;

   };

}
#endif 
