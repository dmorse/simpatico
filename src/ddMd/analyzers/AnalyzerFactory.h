#ifndef DDMD_ANALYZER_FACTORY_H
#define DDMD_ANALYZER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>            // base class template
#include <ddMd/analyzers/Analyzer.h>   // base template parameter

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Factory for DdMd::Analyzer objects.
   *
   * \ingroup DdMd_Factory_Module
   * \ingroup DdMd_Analyzer_Module
   */
   class AnalyzerFactory : public Factory<Analyzer>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation  parent Simulation
      */
      AnalyzerFactory(Simulation& simulation);

      /** 
      * Return pointer to a new Analyzer object.
      *
      * \param className  name of a subclass of Analyzer.
      * \return base class pointer to a new instance of className.
      */
      virtual Analyzer* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent Simulation.
      */
      Simulation& simulation() const
      {  return *simulationPtr_; }

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

   };

}
#endif
