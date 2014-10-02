#ifndef MCMD_LOG_PROGRESS_H
#define MCMD_LOG_PROGRESS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/Analyzer.h>

namespace McMd
{

   using namespace Util;

   /**
   * Periodically write step number to a Log file.
   *
   * \ingroup McMd_Analyzer_Module
   */
   class LogProgress : public Analyzer
   {
   
   public:
  
      /**
      * Constructor.
      */
      LogProgress();

      /**
      * Read interval parameter.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load interval from archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      // using Analyzer::save() and Analyzer::serializable().
   
      /**
      * Write iStep to Log file. 
      *
      * \param iStep MC step index
      */
      void sample(long iStep);
  
   };

}
#endif 
