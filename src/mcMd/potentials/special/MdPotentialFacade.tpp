/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdPotentialFacade.h"
#include <mcMd/simulation/System.h> 

namespace McMd
{

   using namespace Util;

   /**
   * Constructor.
   */
   template <class PotentialType, class FactoryType>
   MdPotentialFacade<PotentialType, FactoryType>::MdPotentialFacade(System& system)
    : MdPotential(false),
      systemPtr_(&system),
      potentialPtr_(0),
      style_()
   {  setClassName("MdPotentialFacade"); }

   /*
   * Destructor.
   */
   template <class PotentialType, class FactoryType>
   MdPotentialFacade<PotentialType, FactoryType>::~MdPotentialFacade()
   {
      if (potentialPtr_) {
         delete potentialPtr_;
      }
   }

   /*
   * Read parameters from file.
   */
   template <class PotentialType, class FactoryType>
   void 
   MdPotentialFacade<PotentialType, FactoryType>::readParameters(std::istream& in)
   {
      read(in, "style", style_);

      FactoryType factory(*systemPtr_);
      potentialPtr_ = factory.factory(style_);
      UTIL_CHECK (potentialPtr_);
      potentialPtr_->readParameters(in);
   }

   /*
   * Compute total eneryg.
   */
   template <class PotentialType, class FactoryType>
   void MdPotentialFacade<PotentialType, FactoryType>::computeEnergy()
   {
      if (!energy_.isSet()) {
         potentialPtr_->computeEnergy();
      }
      energy_.set(potentialPtr_->energy());
   }

   /*
   * Add forces from this potential to all atomic forces.
   */
   template <class PotentialType, class FactoryType>
   void MdPotentialFacade<PotentialType, FactoryType>::addForces()
   {  potentialPtr_->addForces(); }

} 
