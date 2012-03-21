#ifndef DDMD_DIAGNOSTIC_FACTORY_H
#define DDMD_DIAGNOSTIC_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>    // base class template
#include <ddMd/diagnostics/Diagnostic.h>     // base template parameter

namespace DdMd
{

   using namespace Util;

   class System;

   /**
   * Factory for DdMd::Diagnostic objects.
   *
   * \ingroup Factory_Module
   * \ingroup Diagnostic_Module
   */
   class DiagnosticFactory : public Factory<Diagnostic>
   {

   public:

      /**
      * Constructor.
      *
      * \param system     parent System
      */
      DiagnosticFactory(System& system);

      /** 
      * Return pointer to a new Diagnostic object.
      *
      * \param  className name of a subclass of Diagnostic.
      * \return base class pointer to a new instance of className.
      */
      virtual Diagnostic* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent System.
      */
      System& system() const
      {  return *systemPtr_; }

   private:

      /// Pointer to parent System.
      System* systemPtr_;

   };

}
#endif
