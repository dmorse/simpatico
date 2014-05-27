#ifndef MDCF_ANALYZER_FACTORY_H
#define MDCF_ANALYZER_FACTORY_H

/*
* Simpatico - Processor Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>        // base class template
#include <mdCf/analyzers/Analyzer.h>   // base template parameter

namespace MdCf
{

   using namespace Util;

   class Processor;

   /**
   * Factory for MdCf::Analyzer objects.
   *
   * \ingroup MdCf_Analyzer_Module
   */
   class AnalyzerFactory : public Factory<Analyzer>
   {

   public:

      /**
      * Constructor.
      *
      * \param processor     parent Processor
      */
      AnalyzerFactory(Processor& processor);

      /** 
      * Return pointer to a new Analyzer object.
      *
      * \param  className name of a subclass of Analyzer.
      * \return base class pointer to a new instance of className.
      */
      virtual Analyzer* factory(const std::string& className) const;

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
