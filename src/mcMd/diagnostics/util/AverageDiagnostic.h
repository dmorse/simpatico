#ifndef MCMD_AVERAGE_DIAGNOSTIC_H
#define MCMD_AVERAGE_DIAGNOSTIC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <util/accumulators/Average.h>          // member
#include <util/misc/FileMaster.h>  

namespace McMd
{

   using namespace Util;

   /**
   * AverageDiagnostic averages of total potential energy.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   template <class SystemType>
   class AverageDiagnostic : public SystemDiagnostic<SystemType>
   {
   
   public:

      /**   
      * Constructor.
      */
      AverageDiagnostic(SystemType& system);

      /**   
      * Destructor.
      */
      ~AverageDiagnostic();

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
      * Save state to an archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Clear accumulators.
      */
      virtual void setup();

      /**
      * Output results at end of simulation.
      */
      virtual void output();

   protected:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

   };

   /* 
   * Constructor.
   */
   template <class SystemType>
   AverageDiagnostic<SystemType>::AverageDiagnostic(SystemType& system)
    : SystemDiagnostic<SystemType>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1)
   {}

   /* 
   * Constructor.
   */
   template <class SystemType>
   AverageDiagnostic<SystemType>::~AverageDiagnostic()
   {}

   /*
   * Read parameters and initialize.
   */
   template <class SystemType>
   void AverageDiagnostic<SystemType>::readParameters(std::istream& in)
   {

      Diagnostic::readInterval(in);
      Diagnostic::readOutputFileName(in);
      ParamComposite::read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         Diagnostic::fileMaster().openOutputFile(
                             Diagnostic::outputFileName(".dat"), outputFile_);
      }

   }

   /*
   * Load state from an archive.
   */
   template <class SystemType>
   void AverageDiagnostic<SystemType>::loadParameters(
                                       Serializable::IArchive& ar)
   {  
      Diagnostic::loadParameters(ar);
      ParamComposite::loadParameter<int>(ar, "nSamplePerBlock", 
                                         nSamplePerBlock_);
      ar & accumulator_;

      if (nSamplePerBlock_ != accumulator_.nSamplePerBlock()) {
         UTIL_THROW("Inconsistent values of nSamplePerBlock_");
      }

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         Diagnostic::fileMaster().openOutputFile(
                             Diagnostic::outputFileName(".dat"), outputFile_);
      }

   }

   /*
   * Save state to an archive.
   */
   template <class SystemType>
   void AverageDiagnostic<SystemType>::save(Serializable::OArchive& ar)
   {  ar & *this; }


   /*
   * Serialize to/from an archive. 
   */
   template <class SystemType>
   template <class Archive>
   void AverageDiagnostic<SystemType>::serialize(Archive& ar, 
                                                 const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & accumulator_;
   }

   /*
   * Clear accumulators.
   */
   template <class SystemType>
   void AverageDiagnostic<SystemType>::setup()
   {  accumulator_.clear(); }

   /*
   * Output results to file after simulation is completed.
   */
   template <class SystemType>
   void AverageDiagnostic<SystemType>::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }
     
      // Write parameter block to *.prm file
      Diagnostic::fileMaster().openOutputFile(Diagnostic::outputFileName(".prm"), 
                                              outputFile_);
      ParamComposite::writeParam(outputFile_); 
      outputFile_.close();

      // Write average value to *.avefile
      Diagnostic::fileMaster().openOutputFile(Diagnostic::outputFileName(".ave"), 
                                              outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }
   
}
#endif 
