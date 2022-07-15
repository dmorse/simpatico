#ifndef MDPP_ANALYZER_MANAGER_H
#define MDPP_ANALYZER_MANAGER_H

/*
* Simpatico - Processor Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                 // template parameter
#include <util/param/Manager.h>       // base class template

namespace MdPp
{

   using namespace Util;

   /**
   * Manager for a list of Analyzer objects.
   *
   * \ingroup MdPp_Analyzer_Module
   */
   class AnalyzerManager : public Manager<Analyzer>
   {

   public:

      /**
      * Constructor.
      */
      AnalyzerManager();

      /**
      * Destructor.
      */
      virtual ~AnalyzerManager();

      /**
      * Read parameter block (without begin and end).
      *
      * \param in input parameter file stream.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Call setup method of each Analyzer.
      */
      void setup();
 
      /**
      * Call clear method of each Analyzer.
      */
      void clear();
 
      /**
      * Call sample method of each Analyzer, if scheduled.
      *
      * Calls sample methods of each analyzer only if
      * iStep is a multiple of the analyzer interval
      *
      * \param iStep time step counter
      */
      void sample(long iStep);
 
      /**
      * Call output method of each analyzer.
      */
      void output();
 
   };

}
#endif
