#ifndef MCMD_MD_ANALYZER_FACTORY_H
#define MCMD_MD_ANALYZER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>             // base class template
#include <mcMd/analyzers/Analyzer.h>     // base template parameter
#include <mcMd/analyzers/system/SystemAnalyzerFactory.h>  // member

namespace McMd
{

   using namespace Util;

   class MdSimulation;
   class MdSystem;

   /**
   * AnalyzerFactory for an MdSimulation
   *
   * \ingroup McMd_Factory_Module
   * \ingroup McMd_Analyzer_Module
   */
   class MdAnalyzerFactory : public Factory<Analyzer>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      MdAnalyzerFactory(MdSimulation& simulation, MdSystem& system);

      /** 
      * Return pointer to a new Analyzer object.
      *
      * \param  className name of a subclass of Analyzer.
      * \return base class pointer to a new instance of className.
      */
      virtual Analyzer* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent MdSystem.
      */
      MdSystem& system() const
      {  return *systemPtr_; }

      /**
      * Return reference to parent MdSimulation.
      */
      MdSimulation& simulation() const
      {  return *simulationPtr_; }

   private:

      /// Factory for generic System Analyzers.
      SystemAnalyzerFactory systemFactory_;

      /// Pointer to parent MdSimulation.
      MdSimulation* simulationPtr_;

      /// Pointer to parent MdSystem.
      MdSystem* systemPtr_;

   };

}
#endif
