#ifndef MCMD_LOG_PROGRESS_H
#define MCMD_LOG_PROGRESS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/Diagnostic.h>

namespace McMd
{

   using namespace Util;

   /**
   * Periodically write step number to a Log file.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class LogProgress : public Diagnostic
   {
   
   public:
  
      /**
      * Constructor.
      */
      LogProgress();

      /**
      * Read interval.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Write iStep to Log file. 
      *
      * \param iStep MC step index
      */
      void sample(long iStep);
  
   };

}
#endif 
