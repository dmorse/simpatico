#ifndef MDCF_ANALYZER_MANAGER_H
#define MDCF_ANALYZER_MANAGER_H

/*
* Simpatico - Processor Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                 // template parameter
#include <util/param/Manager.h>         // base class template

namespace MdCf
{

   using namespace Util;

   class Processor;

   /**
   * Manager for a list of Analyzer objects.
   *
   * \ingroup MdCf_Manager_Module
   * \ingroup MdCf_Analyzer_Module
   */
   class AnalyzerManager : public Manager<Analyzer>
   {

   public:

      /**
      * Constructor.
      */
      AnalyzerManager(Processor& processor);

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

      /**
      * Return pointer to a new default factory.
      */
      virtual Factory<Analyzer>* newDefaultFactory() const;

   private:

      /// Pointer to parent Processor.
      Processor* processorPtr_;
 
   };

}
#endif
