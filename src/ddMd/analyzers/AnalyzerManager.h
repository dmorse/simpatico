#ifndef DDMD_ANALYZER_MANAGER_H
#define DDMD_ANALYZER_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                 // template parameter
#include <util/param/Manager.h>         // base class template

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Manager for a list of Analyzer objects.
   *
   * \ingroup DdMd_Manager_Module
   * \ingroup DdMd_Analyzer_Module
   */
   class AnalyzerManager : public Manager<Analyzer>
   {

   public:

      /**
      * Constructor.
      */
      AnalyzerManager(Simulation& simulation);

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
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
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
      * This function should be called every time step during
      * a simulation. It executes a loop over analyzers if the
      * static member Analyzer::baseInterval is positive on 
      * time steps for which iStep is an interval multiple of 
      * Analyzer::baseInterval. It call the sample() function
      * of each Analyzer object only when iStep is an integer
      * multiple of the interval variable for that Analyzer.
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

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;
 
   };

}
#endif
