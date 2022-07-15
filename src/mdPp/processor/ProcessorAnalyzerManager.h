#ifndef MDPP_PROCESSOR_ANALYZER_MANAGER_H
#define MDPP_PROCESSOR_ANALYZER_MANAGER_H

/*
* Simpatico - Processor Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mdPp/analyzers/AnalyzerManager.h> // base class

namespace MdPp
{

   using namespace Util;

   class Processor;

   /**
   * Manager for a list of Analyzer objects.
   *
   * \ingroup MdPp_Analyzer_Module
   */
   class ProcessorAnalyzerManager : public AnalyzerManager
   {

   public:

      /**
      * Constructor.
      */
      ProcessorAnalyzerManager(Processor& processor);

      /**
      * Destructor.
      */
      virtual ~ProcessorAnalyzerManager();

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
