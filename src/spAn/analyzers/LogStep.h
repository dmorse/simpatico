#ifndef SPAN_LOG_STEP_H
#define SPAN_LOG_STEP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/analyzers/Analyzer.h>       // base class template

namespace SpAn
{

   using namespace Util;

   /**
   * Write step number to log.
   *
   * \ingroup SpAn_Analyzer_Module
   */
   class LogStep : public Analyzer
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param processor reference to parent Processor
      */
      LogStep(Processor &processor)
       : Analyzer(processor)
      {}
  
      /** 
      * Read parameters from file.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in){}
  
      /** 
      * Determine number of molecules and allocate memory.
      */
      virtual void setup(){};
   
      /** 
      * Evaluate end-to-end vectors of all chains, add to ensemble.
      *
      * \param iStep counter for number of steps
      */
      virtual void sample(long iStep)
      {  std::cout << "iStep = " << iStep << std::endl; }

      /**
      * Output results to file after simulation is completed.
      */
      virtual void output(){};
   
   };

}
#endif
