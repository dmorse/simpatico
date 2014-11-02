#ifndef MCMD_SYSTEM_ANALYZER_H
#define MCMD_SYSTEM_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/Analyzer.h>
#include <mcMd/simulation/Simulation.h>
#include <util/global.h>                   // assert in inline function

namespace McMd
{

   using namespace Util;

   /**
   * Template for Analyzer associated with one System.
   *
   * The FileMaster associated with a SystemAnalyzer is the one
   * used by the parent System.
   *
   * \ingroup McMd_Analyzer_Module
   */
   template <class SystemType>
   class SystemAnalyzer : public Analyzer 
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent System
      */
      SystemAnalyzer(SystemType& system);
  
      /**
      * Destructor.
      */
      virtual ~SystemAnalyzer();
  
   protected:
  
      /** 
      * Return reference to parent system.
      */
      SystemType& system();

   private:

      /// Pointer to parent system.
      SystemType* systemPtr_;

   };

   /*
   * Constructor.
   */
   template <class SystemType>
   SystemAnalyzer<SystemType>::SystemAnalyzer(SystemType& system) 
    : Analyzer(),
      systemPtr_(&system)
   {  setFileMaster(system.fileMaster()); }

   /*
   * Destructor.
   */
   template <class SystemType>
   SystemAnalyzer<SystemType>::~SystemAnalyzer() 
   {}

   /*
   * Return reference to parent system.
   */
   template <class SystemType>
   inline SystemType& SystemAnalyzer<SystemType>::system()
   {
      assert(systemPtr_);
      return *systemPtr_; 
   }

}
#endif
