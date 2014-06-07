#ifndef DDMD_SP_ANALYZER_MANAGER_H
#define DDMD_SP_ANALYZER_MANAGER_H

/*
* Simpatico - Processor Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpAnalyzer.h"                 // template parameter
#include <util/param/Manager.h>         // base class template

namespace DdMd
{

   using namespace Util;

   class Processor;

   /**
   * Manager for a list of SpAnalyzer objects.
   *
   * \ingroup DdMd_Sp_Analyzer_Module
   */
   class SpAnalyzerManager : public Manager<SpAnalyzer>
   {

   public:

      /**
      * Constructor.
      */
      SpAnalyzerManager(Processor& processor);

      /**
      * Destructor.
      */
      virtual ~SpAnalyzerManager();

      /**
      * Read parameter block (without begin and end).
      *
      * \param in input parameter file stream.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Call setup method of each SpAnalyzer.
      */
      void setup();
 
      /**
      * Call clear method of each SpAnalyzer.
      */
      void clear();
 
      /**
      * Call sample method of each SpAnalyzer, if scheduled.
      *
      * Calls sample methods of each analyzer only if
      * iStep is a multiple of the analyzer interval
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
      virtual Factory<SpAnalyzer>* newDefaultFactory() const;

   private:

      /// Pointer to parent Processor.
      Processor* processorPtr_;
 
   };

}
#endif
