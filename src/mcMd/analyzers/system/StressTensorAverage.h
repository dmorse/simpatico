#ifndef MCMD_STRESS_TENSOR_AVERAGE_H
#define MCMD_STRESS_TENSOR_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/space/Tensor.h>
#include <util/accumulators/Average.h>     // member template

namespace McMd
{

   using namespace Util;

   /**
   * Periodically write (tensor) StressTensor to file.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   template <class SystemType>
   class StressTensorAverage : public SystemAnalyzer<SystemType>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent SystemType object. 
      */
      StressTensorAverage(SystemType& system);
   
      /**
      * Destructor.
      */
      virtual ~StressTensorAverage()
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
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);
  
      /**
      * Clear nSample counter.
      */
      virtual void clear();

      /**
      * Setup accumulator!
      */
      virtual void setup();
  
      /**
      * Sample virial stress components to accumulators
      *
      * \param iStep MD or MC step index
      */
      virtual void sample(long iStep);

      /**
      * Output final results to file
      */
      virtual void output();

   private:
 
      /// Output file stream
      std::ofstream  outputFile_;
      
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

      /// Has readParam been called?
      long  isInitialized_;

      using SystemAnalyzer<SystemType>::readInterval;
      using SystemAnalyzer<SystemType>::readOutputFileName;
      using SystemAnalyzer<SystemType>::read;
      using SystemAnalyzer<SystemType>::writeParam;
      using SystemAnalyzer<SystemType>::loadParameter;
      using SystemAnalyzer<SystemType>::isAtInterval;
      using SystemAnalyzer<SystemType>::outputFileName;
      using SystemAnalyzer<SystemType>::fileMaster;
      using SystemAnalyzer<SystemType>::system;
   
   };

   /*
   * Constructor.
   */
   template <class SystemType>
   StressTensorAverage<SystemType>::StressTensorAverage(SystemType& system)
    : SystemAnalyzer<SystemType>(system),
      outputFile_(),
      sxxAccumulator_(),
      sxyAccumulator_(),
      sxzAccumulator_(),
      syxAccumulator_(),
      syyAccumulator_(),
      syzAccumulator_(),
      szxAccumulator_(),
      szyAccumulator_(),
      szzAccumulator_(),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {}

   /*
   * Read parameters and initialize.
   */
   template <class SystemType>
   void StressTensorAverage<SystemType>::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      sxxAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      sxyAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      sxzAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      syxAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      syyAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      syzAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      szxAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      szyAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      szzAccumulator_.setNSamplePerBlock(nSamplePerBlock_);

      std::string filename;
      filename  = outputFileName();
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   template <class SystemType>
   void StressTensorAverage<SystemType>::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);  
      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);

      sxxAccumulator_.loadParameters(ar);
      sxyAccumulator_.loadParameters(ar);
      sxzAccumulator_.loadParameters(ar);
      syxAccumulator_.loadParameters(ar);
      syyAccumulator_.loadParameters(ar);
      syzAccumulator_.loadParameters(ar);
      szxAccumulator_.loadParameters(ar);
      szyAccumulator_.loadParameters(ar);
      szzAccumulator_.loadParameters(ar);

      if (nSamplePerBlock_) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   template <class SystemType>
   void StressTensorAverage<SystemType>::save(Serializable::OArchive& ar)
   { ar & *this; }


   /*
   * Serialize to/from an archive. 
   */
   template <class SystemType>
   void StressTensorAverage<SystemType>::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }

      sxxAccumulator_.clear();
      sxyAccumulator_.clear();
      sxzAccumulator_.clear();
      syxAccumulator_.clear();
      syyAccumulator_.clear();
      syzAccumulator_.clear();
      szxAccumulator_.clear();
      szyAccumulator_.clear();
      szzAccumulator_.clear();       
   }

   /* 
   * Evaluate pressure, and add to accumulator.
   */
   template <class SystemType>
   void StressTensorAverage<SystemType>::sample(long iStep)
   {
      if (isAtInterval(iStep)){

         SystemType& sys=system(); 

         Tensor total;
         Tensor virial;
         Tensor kinetic;

         sys.computeVirialStress(total);
         sys.computeKineticStress(kinetic);
         total.add(virial, kinetic);

         sxxAccumulator_.sample(total(0,0));
         sxyAccumulator_.sample(total(0,1));
         sxzAccumulator_.sample(total(0,2));
         syxAccumulator_.sample(total(1,0));
         syyAccumulator_.sample(total(1,1));
         syzAccumulator_.sample(total(1,2));
         szxAccumulator_.sample(total(2,0));
         szyAccumulator_.sample(total(2,1));
         szzAccumulator_.sample(total(2,2));
     }
   }

   /*
   * Output results to file after simulation is completed.
   */
   template <class SystemType>
   void StressTensorAverage<SystemType>::output() 
   {
      outputFile_ << 
                  "Sxx=" << Dbl(sxxAccumulator_.average(), 17)<< "  +-  " << Dbl(sxxAccumulator_.error(), 9, 8) << "\n" << 
                  "Sxy=" << Dbl(sxyAccumulator_.average(), 17)<< "  +-  " << Dbl(sxyAccumulator_.error(), 9, 8) << "\n" << 
                  "Sxz=" << Dbl(sxzAccumulator_.average(), 17)<< "  +-  " << Dbl(sxzAccumulator_.error(), 9, 8) << "\n" << 
                  "Syx=" << Dbl(syxAccumulator_.average(), 17)<< "  +-  " << Dbl(syxAccumulator_.error(), 9, 8) << "\n" << 
                  "Syy=" << Dbl(syyAccumulator_.average(), 17)<< "  +-  " << Dbl(syyAccumulator_.error(), 9, 8) << "\n" << 
                  "Syz=" << Dbl(syzAccumulator_.average(), 17)<< "  +-  " << Dbl(syzAccumulator_.error(), 9, 8) << "\n" << 
                  "Szx=" << Dbl(szxAccumulator_.average(), 17)<< "  +-  " << Dbl(szxAccumulator_.error(), 9, 8) << "\n" << 
                  "Szy=" << Dbl(szyAccumulator_.average(), 17)<< "  +-  " << Dbl(szyAccumulator_.error(), 9, 8) << "\n" << 
                  "Szz=" << Dbl(szzAccumulator_.average(), 17)<< "  +-  " << Dbl(szzAccumulator_.error(), 9, 8) << "\n" << 
                  std::endl;

      outputFile_.close(); 
   }

}
#endif 
