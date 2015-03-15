#ifndef MCMD_MC_ANALYZER_MANAGER_H
#define MCMD_MC_ANALYZER_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/AnalyzerManager.h>  // base class

namespace McMd
{

   using namespace Util;

   class Simulation;
   class McSimulation;
   class McSystem;

   /**
   * Manager for Analyzer objects in an McSimulation.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_Analyzer_Mc_Module
   */
   class McAnalyzerManager : public AnalyzerManager
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      */
      McAnalyzerManager(McSimulation& simulation);

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      * \param system     associated McSystem
      */
      McAnalyzerManager(McSimulation& simulation, McSystem& system);

   protected:

      /**
      * Instantiate and return pointer to a new McAnalyzerFactory object.
      */
      virtual Factory<Analyzer>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation.
      McSimulation* simulationPtr_;

      /// Pointer to associated McSystem.
      McSystem*   systemPtr_;

   };

}
#endif
