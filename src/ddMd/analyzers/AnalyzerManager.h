#ifndef DDMD_ANALYZER_MANAGER_H
#define DDMD_ANALYZER_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
      * Calls sample methods of each analyzer only if:
      * - Analyzer::baseInterval is positive
      * - iStep is a multiple of Analyzer::baseInterval
      * - iStep is a multiple of the analyzer interval
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
