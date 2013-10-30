#ifndef MCMD_MC_ANALYZER_FACTORY_H
#define MCMD_MC_ANALYZER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>     // base class template
#include <mcMd/analyzers/Analyzer.h>    // base class template parameter
#include <mcMd/analyzers/system/SystemAnalyzerFactory.h>    // member

namespace McMd
{

   using namespace Util;

   class McSimulation;
   class McSystem;

   /**
   * AnalyzerFactory for an McSimulation
   *
   * \ingroup McMd_Factory_Module
   * \ingroup McMd_Analyzer_Module
   */
   class McAnalyzerFactory : public Factory<Analyzer>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      * \param system     parent McSystem
      */
      McAnalyzerFactory(McSimulation& simulation, McSystem& system);

      /** 
      * Return pointer to a new Analyzer object.
      *
      * \param  className name of a subclass of Analyzer.
      * \return base class pointer to a new instance of className.
      */
      virtual Analyzer* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent McSystem.
      */
      McSystem& system() const
      { return *systemPtr_; }

      /**
      * Return reference to parent McSimulation.
      */
      McSimulation& simulation() const
      { return *simulationPtr_; }

   private:

      // Factory for analyzers for any System
      SystemAnalyzerFactory systemFactory_;

      // Pointer to parent Simulation
      McSimulation *simulationPtr_;

      // Pointer to parent System
      McSystem  *systemPtr_;

   };

}
#endif
