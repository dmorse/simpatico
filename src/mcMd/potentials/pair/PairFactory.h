#ifndef INTER_NOPAIR
#ifndef MCMD_PAIR_FACTORY_H
#define MCMD_PAIR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <iostream>
#include <vector>

namespace McMd
{

   class System;
   class MdPairPotential;
   class McPairPotential;

   /**
   * Factory for subclasses MdPairPotential or McPairPotential.
   *
   * This class is not derived from the Factory template because it
   * provides a different interface, with separate mcFactory() and
   * mdFactory() methods to create McPairPotential and MdPairPotential
   * objects. Unlike a standard Factory, it also does not provide a 
   * readObject method, and so does not require a param communicator.
   * 
   * \ingroup McMd_Pair_Module
   */
   class PairFactory 
   {
   
   public:
   
      /**
      * Constructor.
      */
      PairFactory();

      /**
      * Destructor.
      */
      virtual ~PairFactory()
      {}

      /**
      * Add a new subfactory to the list.
      *
      * If this Factory has a param communicator, this method also sets 
      * the child subfactory communicator to that of the parent.
      *
      * \param subfactory New subfactory to be added
      */
      void addSubfactory(PairFactory& subfactory);

      /**
      * Return a pointer to a new McPairPotential, if possible.
      *
      * \param subclass name of desired subclass of McPairPotential
      * \param system associated System
      */
      virtual McPairPotential* mcFactory(const std::string& subclass, System& system) const;

      /**
      * Return a pointer to a new McPairPotential, if possible.
      *
      * \param subclass name of desired subclass of MdPairPotential
      * \param system   associated System
      */
      virtual MdPairPotential* mdFactory(const std::string& subclass, System& system) const;

      /**
      * Create an MdPairPotential from a McPairPotential.
      *
      * \param potential McPairPotential to be cloned
      */
      virtual MdPairPotential* mdFactory(McPairPotential& potential) const;

   protected:

      /**
      * Search subfactories for match to McPairPotential subclass name. 
      *
      * \param className name of subclass
      * \param system    associated System
      * \return base class pointer to new McPairPotential, or a null pointer.
      */
      McPairPotential* tryMcSubfactories(const std::string& className, System& system) const;

      /**
      * Search subfactories for match to MdPairPotential subclass name. 
      *
      * \param  className name of subclass
      * \param  system    associated System
      * \return base class pointer to new MdPairPotential, or a null pointer.
      */
      MdPairPotential* tryMdSubfactories(const std::string& className, System& system) const;

      /**
      * Search through subfactories for match for an mc potential.
      *
      * \param  potential McPairPotential to be cloned
      * \return base class pointer to new MdPairPotential, or a null pointer.
      */
      MdPairPotential* tryMdSubfactories(McPairPotential& potential) const;

   private:

      /// Vector of pointers to child PairFactory objects.
      std::vector< PairFactory* > subfactories_;

   };
 
}
#endif
#endif
