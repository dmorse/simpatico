#ifndef DDMD_DIAGNOSTIC_MANAGER_H
#define DDMD_DIAGNOSTIC_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Diagnostic.h"                 // template parameter
#include <util/param/Manager.h>         // base class template

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Manager for a list of Diagnostic objects.
   *
   * \ingroup DdMd_Manager_Module
   * \ingroup DdMd_Diagnostic_Module
   */
   class DiagnosticManager : public Manager<Diagnostic>
   {

   public:

      /**
      * Constructor.
      */
      DiagnosticManager(Simulation& simulation);

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
      * Call setup method of each Diagnostic.
      */
      void setup();
 
      /**
      * Call clear method of each Diagnostic.
      */
      void clear();
 
      /**
      * Call sample method of each Diagnostic, if scheduled.
      *
      * Calls sample methods of each diagnostic only if:
      * - Diagnostic::baseInterval is positive
      * - iStep is a multiple of Diagnostic::baseInterval
      * - iStep is a multiple of the diagnostic interval
      *
      * \param iStep time step counter
      */
      void sample(long iStep);
 
      /**
      * Call output method of each diagnostic.
      */
      void output();

      /**
      * Return pointer to a new default factory.
      */
      virtual Factory<Diagnostic>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;
 
   };

}
#endif
