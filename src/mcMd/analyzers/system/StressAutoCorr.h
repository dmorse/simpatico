#ifndef MCMD_STRESS_AUTO_CORR_H
#define MCMD_STRESS_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>
#include <util/space/Tensor.h>
#include <util/accumulators/AutoCorrelation.tpp>     // member template

namespace McMd
{

   using namespace Util;

   /**
   * Compute stress autocorrelation for an isotropic system
   *
   * Typename SystemType can be McSystem or MdSystem.
   *
   * \ingroup McMd_Analyzer_Module
   */
   template <class SystemType>
   class StressAutoCorr : public SystemAnalyzer<SystemType>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent SystemType object. 
      */
      StressAutoCorr(SystemType& system);
   
      /**
      * Destructor.
      */
      virtual ~StressAutoCorr();
   
      /**
      * Read parameters from file and initialize.
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
      * Setup before beginning of loop.
      */
      virtual void setup();
  
      /**
      * Sample stress tensor and add to accumulator.
      *
      * \param iStep MD or MC step index
      */
      virtual void sample(long iStep);

      /**
      * Output results to file (parameters and autocorrelation function).
      */
      virtual void output();

   protected:

      /**
      * Compute total stress tensor.
      * 
      * \param stress Stress tensor (on return).
      */ 
      virtual void computeStress(Tensor& stress) = 0;

      using SystemAnalyzer<SystemType>::readInterval;
      using SystemAnalyzer<SystemType>::readOutputFileName;
      using SystemAnalyzer<SystemType>::read;
      using SystemAnalyzer<SystemType>::writeParam;
      using SystemAnalyzer<SystemType>::loadParameter;
      using SystemAnalyzer<SystemType>::isAtInterval;
      using SystemAnalyzer<SystemType>::outputFileName;
      using SystemAnalyzer<SystemType>::fileMaster;
      using SystemAnalyzer<SystemType>::system;
  
   private:
 
      /// Output file stream
      std::ofstream  outputFile_;
 
      /// Statistical accumulator.
      AutoCorrelation<Tensor, double>  accumulator_;

      /// Number of samples per block average output
      int  capacity_;

      /// Maximum id for AutoCorrStage
      int  maxStageId_;

      /// BlockFactor for AutoCorrStage algorithm
      int  blockFactor_;

      /// Has readParam been called?
      long  isInitialized_;

   };

   /*
   * Constructor.
   */
   template <class SystemType>
   StressAutoCorr<SystemType>::StressAutoCorr(SystemType& system)
    : SystemAnalyzer<SystemType>(system),
      outputFile_(),
      accumulator_(),
      capacity_(64),
      maxStageId_(10),
      blockFactor_(2),
      isInitialized_(false)
   {}

   /*
   * Destructor.
   */
   template <class SystemType>
   StressAutoCorr<SystemType>::~StressAutoCorr()
   {} 
   
   /*
   * Read parameters and initialize.
   */
   template <class SystemType>
   void StressAutoCorr<SystemType>::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read(in, "capacity", capacity_);
      bool isRequired = false;
      read(in, "maxStageId", maxStageId_, isRequired);

      accumulator_.setParam(capacity_, maxStageId_, blockFactor_);
      accumulator_.clear();

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   template <class SystemType>
   void StressAutoCorr<SystemType>::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);

      loadParameter(ar, "capacity", capacity_);
      ar & accumulator_;

      if (accumulator_.bufferCapacity() != capacity_) {
         UTIL_THROW("Inconsistent values of capacity");
      }

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   template <class SystemType>
   void StressAutoCorr<SystemType>::save(Serializable::OArchive& ar)
   {  ar & *this; }


   /*
   * Serialize to/from an archive. 
   */
   template <class SystemType>
   template <class Archive>
   void 
   StressAutoCorr<SystemType>::serialize(Archive& ar, 
                                      const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & capacity_;
      ar & accumulator_;
   }

   /*
   * Set up immediately before simulation.
   */
   template <class SystemType>
   void StressAutoCorr<SystemType>::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized");
      }  
      accumulator_.clear(); 
   }

   /* 
   * Evaluate pressure, and add to accumulator.
   */
   template <class SystemType>
   void StressAutoCorr<SystemType>::sample(long iStep)
   {
      if (isAtInterval(iStep)){

         #if 0
         // Compute stress
         Tensor total;
         Tensor virial;
         Tensor kinetic;
         system().computeVirialStress(virial);
         system().computeKineticStress(kinetic);
         total.add(virial, kinetic);
         #endif

         Tensor stress;
         computeStress(stress);

         // Remove trace
         double pressure = 0.0;
         int i, j;
         for (i = 0; i < Dimension; ++i) {
            pressure += stress(i,i);
         }
         pressure = pressure/double(Dimension);
         for (i = 0; i < Dimension; ++i) {
            stress(i,i) -= pressure;
         }

         double factor = 1.0/sqrt(10.0);
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < Dimension; ++j) {
               stress(i,j) *= factor;
            }
         }

         accumulator_.sample(stress);
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   template <class SystemType>
   void StressAutoCorr<SystemType>::output() 
   {
      // Write parameters to *.prm file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_ << std::endl;
      outputFile_ << "bufferCapacity  " 
                  << accumulator_.bufferCapacity() << std::endl;
      outputFile_ << "nSample         " 
                  << accumulator_.nSample() << std::endl;
      outputFile_ << std::endl;
      outputFile_.close();

      // Write autocorrelation function to data file
      fileMaster().openOutputFile(outputFileName(".corr"), outputFile_);
      accumulator_.output(outputFile_);
      outputFile_.close(); 
   }

}
#endif 
