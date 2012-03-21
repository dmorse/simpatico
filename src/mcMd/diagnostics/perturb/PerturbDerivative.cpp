#ifdef MCMD_PERTURB
#ifndef PERTURB_DERIVATIVE_CPP
#define MCMD_PERTURB_DERIVATIVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbDerivative.h"           // class header
#include <mcMd/perturb/Perturbation.h>  
#include <mcMd/util/FileMaster.h>  

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   PerturbDerivative::PerturbDerivative(System& system)
    : SystemDiagnostic<System>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1),
      parameterIndex_(0)
   {}

   /*
   * Read parameters and initialize.
   */
   void PerturbDerivative::readParam(std::istream& in)
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      read<int>(in,"parameterIndex", parameterIndex_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

   }

   /* 
   * Evaluate perturb derivative, and add to accumulator.
   */
   void PerturbDerivative::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         accumulator_.sample(system().perturbation().derivative(parameterIndex_), outputFile_);
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void PerturbDerivative::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }
   
}
#endif    // ifndef PERTURB_DERIVATIVE_CPP
#endif    // ifdef  MCMD_PERTURB
