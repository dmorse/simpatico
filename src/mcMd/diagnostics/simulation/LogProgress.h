#ifndef MCMD_LOG_PROGRESS_H
#define MCMD_LOG_PROGRESS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/Diagnostic.h>

namespace McMd
{

   using namespace Util;

   /**
   * Periodically write step number to a Log file.
   *
   * \ingroup Diagnostic_Module
   */
   class LogProgress : public Diagnostic
   {
   
   public:
   
      /**
      * Read interval.
      *
      * \param in input parameter file
      */
      virtual void readParam(std::istream& in);
   
      /**
      * Write iStep to Log file. 
      *
      * \param iStep MC step index
      */
      void sample(long iStep);
  
   };

}
#endif 
