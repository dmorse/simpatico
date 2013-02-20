#ifndef MCMD_DIAGNOSTIC_MANAGER_H
#define MCMD_DIAGNOSTIC_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Diagnostic.h"                  // template parameter
#include <util/param/Manager.h>          // base class template

namespace McMd
{

   using namespace Util;

   /**
   * Manager for a list of Diagnostic objects.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_Diagnostic_Module
   */
   class DiagnosticManager : public Manager<Diagnostic>
   {

   public:

      /**
      * Constructor.
      */
      DiagnosticManager();

      /**
      * Destructor.
      */
      virtual ~DiagnosticManager();

      /**
      * Read parameter file. 
      *
      * \param in input parameter file stream.
      */
      virtual void readParam(std::istream &in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Call initialize method of each Diagnostic.
      * 
      * This method should be called just before the main
      * simulation loop, after an initial configuration is
      * known.
      */
      void setup();
 
      /**
      * Call sample method of each Diagnostic.
      */
      void sample(long iStep);
 
      /**
      * Call output method of each diagnostic.
      */
      void output();

   };

}
#endif
