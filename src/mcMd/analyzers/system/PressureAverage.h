#ifndef MCMD_PRESSURE_AVERAGE_H
#define MCMD_PRESSURE_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <util/accumulators/Average.h>          // member
#include <util/misc/FileMaster.h>  
#include <util/archives/Serializable_includes.h>

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * PressureAverage evaluates average pressure.
   *
   * The SystemType may be McSystem and MdSystem. The use of a template
   * is possible because both types of system use a similar interface.
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   template <class SystemType>
   class PressureAverage : public SystemAnalyzer<SystemType>
   {
   
   public:

      /**   
      * Constructor.
      */
      PressureAverage(SystemType& system);

      /**   
      * Destructor.
      */
      ~PressureAverage();

      /**
      * Read parameters and initialize.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Set up before simulation.
      */
      virtual void setup();
   
      /* 
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Has readParam been called?
      bool    isInitialized_;

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
   PressureAverage<SystemType>::PressureAverage(SystemType& system)
    : SystemAnalyzer<SystemType>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(-1),
      isInitialized_(false)
   {}

   /*
   * Destructor.
   */
   template <class SystemType>
   PressureAverage<SystemType>::~PressureAverage()
   {}

   /*
   * Read parameters and initialize.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // Open output file for block averages, if nSamplePerBlock != 0.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      loadParameter(ar, "nSamplePerBlock", nSamplePerBlock_);
      ar & accumulator_;

      if (accumulator_.nSamplePerBlock() != nSamplePerBlock_) {
         UTIL_THROW("Inconsistent values of nSamplePerBlock");
      }

      // Open output file for block averages, if nSamplePerBlock != 0.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }


   /*
   * Save state to archive.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Serialize to/from an archive. 
   */
   template <class SystemType>
   template <class Archive>
   void PressureAverage<SystemType>::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & accumulator_;
   }

   /*
   * Set up immediately before simulation.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::setup() 
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
   void PressureAverage<SystemType>::sample(long iStep) 
   {
      double pressure;
      SystemType* systemPtr_ = &system();
      systemPtr_->computeStress(pressure);
      if (isAtInterval(iStep)) {
         accumulator_.sample(pressure, outputFile_);
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::output() 
   {
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      // Write parameters to *.prm file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Write average and variance to *.ave file
      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
#endif
