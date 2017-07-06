#ifndef MCMD_SPECIAL_FACTORY_H
#define MCMD_SPECIAL_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <iostream>
#include <vector>

namespace McMd
{

   class System;
   class MdPotential;

   /**
   * Factory for specialized subclasses of MdPotential.
   *
   * \ingroup McMd_Special_Module
   */
   class SpecialFactory 
   {
   
   public:
   
      /**
      * Constructor.
      */
      SpecialFactory();

      /**
      * Destructor.
      */
      virtual ~SpecialFactory()
      {}

      /**
      * Add a new subfactory to the list.
      *
      * If this Factory has a param communicator, this method also sets 
      * the child subfactory communicator to that of the parent.
      *
      * \param subfactory New subfactory to be added
      */
      void addSubfactory(SpecialFactory& subfactory);

      /**
      * Return a pointer to a new MdPotential, if possible.
      *
      * \param subclass name of desired subclass of MdPotential
      * \param system   associated System
      */
      virtual MdPotential* mdFactory(const std::string& subclass, System& system) const;

   protected:

      /**
      * Search subfactories for match to MdPotential subclass name. 
      *
      * \param  className name of subclass
      * \param  system    associated System
      * \return base class pointer to new MdPotential, or a null pointer.
      */
      MdPotential* tryMdSubfactories(const std::string& className, System& system) const;

   private:

      /// Vector of pointers to child SpecialFactory objects.
      std::vector< SpecialFactory* > subfactories_;

   };
 
}
#endif

