#ifndef DDMD_SP_ANALYZER_FACTORY_H
#define DDMD_SP_ANALYZER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>        // base class template
#include <ddMd/sp/analyzers/SpAnalyzer.h>   // base template parameter

namespace DdMd
{

   using namespace Util;

   class Processor;

   /**
   * Factory for DdMd::SpAnalyzer objects.
   *
   * \ingroup DdMd_Sp_Analyzer_Module
   */
   class SpAnalyzerFactory : public Factory<SpAnalyzer>
   {

   public:

      /**
      * Constructor.
      *
      * \param processor associated physical Processor
      */
      SpAnalyzerFactory(Processor& processor);

      /** 
      * Return pointer to a new SpAnalyzer object.
      *
      * \param  className name of a subclass of SpAnalyzer.
      * \return base class pointer to a new instance of className.
      */
      virtual SpAnalyzer* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent Processor.
      */
      Processor& processor() const
      {  return *processorPtr_; }

   private:

      /// Pointer to parent Processor.
      Processor* processorPtr_;

   };

}
#endif
