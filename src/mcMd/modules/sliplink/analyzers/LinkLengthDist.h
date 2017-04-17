#ifndef LINK_LENGTH_DIST_H
#define LINK_LENGTH_DIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>    // base class template
#include <mcMd/simulation/System.h>               // base class template parameter
#include <util/accumulators/Distribution.h>
#include <util/containers/DArray.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * LinkLengthDist evaluates the distribution function of the lengths of the links.
   *
   * \ingroup Analyzer_Module
   */
   class LinkLengthDist : public SystemAnalyzer<System>
   {
   
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      LinkLengthDist(System &system);

      /**
      * Read parameters from file.  
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /** 
      * Add particle pairs to LinkLengthDist histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /** 
      * Output results to output file.
      */
      virtual void output();

   private:

      // Output file stream
      std::ofstream outputFile_;

      // Distribution statistical accumulator
      Distribution  accumulator_;

      /// Index of relevant Species.
      int     speciesId_;
  

   };

}
#endif
