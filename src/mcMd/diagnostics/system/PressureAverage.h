#ifndef MCMD_PRESSURE_AVERAGE_H
#define MCMD_PRESSURE_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <util/accumulators/Average.h>          // member
#include <mcMd/util/FileMaster.h>  
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
   * \ingroup Diagnostic_Module
   */
   template <class SystemType>
   class PressureAverage : public SystemDiagnostic<SystemType>
   {
   
   public:

      /**   
      * Constructor.
      */
      PressureAverage(SystemType& system);

      /**
      * Read parameters and initialize.
      */
      virtual void readParam(std::istream& in);

      /** 
      * Set up before simulation.
      */
      virtual void initialize();
   
      /* 
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

      /**
      * Save state to binary file archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchiveType& ar);

      /**
      * Load state from a binary file archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void load(Serializable::IArchiveType& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Has readParam been called?
      bool    isInitialized_;

      using SystemDiagnostic<SystemType>::readInterval;
      using SystemDiagnostic<SystemType>::readOutputFileName;
      using SystemDiagnostic<SystemType>::writeParam;
      using SystemDiagnostic<SystemType>::isAtInterval;
      using SystemDiagnostic<SystemType>::outputFileName;
      using SystemDiagnostic<SystemType>::fileMaster;
      using SystemDiagnostic<SystemType>::system;

   };

   /*
   * Constructor.
   */
   template <class SystemType>
   PressureAverage<SystemType>::PressureAverage(SystemType& system)
    : SystemDiagnostic<SystemType>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(-1),
      isInitialized_(false)
   {}

   /*
   * Read parameters and initialize.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::readParam(std::istream& in)
   {

      readInterval(in);
      readOutputFileName(in);
      SystemDiagnostic<SystemType>::read(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      isInitialized_ = true;
   }

   /*
   * Set up immediately before simulation.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::initialize() 
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

   /*
   * Save state to binary file archive.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   template <class SystemType>
   void PressureAverage<SystemType>::load(Serializable::IArchiveType& ar)
   { ar & *this; }

   /*
   * Serialize to/from an archive. 
   */
   template <class SystemType>
   template <class Archive>
   void PressureAverage<SystemType>::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }
      ar & accumulator_;
   }

}
#endif
