#ifndef AVERAGE_DIAGNOSTIC_H
#define AVERAGE_DIAGNOSTIC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <util/accumulators/Average.h>          // member
#include <mcMd/util/FileMaster.h>  

namespace McMd
{

   using namespace Util;

   /**
   * AverageDiagnostic averages of total potential energy.
   *
   * \ingroup Diagnostic_Module
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
      * Read parameters and initialize.
      */
      virtual void readParam(std::istream& in);

      /**
      * Clear accumulators.
      */
      virtual void initialize();

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
   * Read parameters and initialize.
   */
   template <class SystemType>
   void AverageDiagnostic<SystemType>::readParam(std::istream& in)
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
   * Clear accumulators.
   */
   template <class SystemType>
   void AverageDiagnostic<SystemType>::initialize()
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
