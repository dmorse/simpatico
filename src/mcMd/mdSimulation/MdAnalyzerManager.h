#ifndef MCMD_MD_ANALYZER_MANAGER_H
#define MCMD_MD_ANALYZER_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/AnalyzerManager.h>  // base class

namespace Util { template <typename T> class Factory; }

namespace McMd
{

   using namespace Util;

   class Simulation;
   class MdSimulation;
   class MdSystem;

   /**
   * Manager for Analyzer objects in an MdSimulation.
   *
   * The default Factory<Analyzer> object for an MdAnalyzerManager is
   * an MdAnalyzerFactory.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_Analyzer_Module
   */
   class MdAnalyzerManager : public AnalyzerManager
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent MdSimulation
      */
      MdAnalyzerManager(MdSimulation& simulation);

      /**
      * Constructor.
      *
      * \param simulation parent Simulation
      * \param system     associated MdSystem
      */
      MdAnalyzerManager(MdSimulation& simulation, MdSystem& system);

      /**
      * Destructor.
      */
      virtual ~MdAnalyzerManager();

   protected:

      /**
      * Create and return pointer to a new MdAnalyzerFactory object
      */
      virtual Util::Factory<Analyzer>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation
      MdSimulation* simulationPtr_;

      /// Pointer to associated MdSystem.
      MdSystem*   systemPtr_;

   };

}
#endif
