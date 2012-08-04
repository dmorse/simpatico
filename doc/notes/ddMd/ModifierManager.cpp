#ifndef DDMD_MODIFIER_MANAGER_CPP
#define DDMD_MODIFIER_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ModifierManager.h" 
//#include "ModifierFactory.h" 

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ModifierManager::ModifierManager(Simulation& simulation)
   : Manager<Modifier>(),
     simulationPtr_(&simulation)
   {}

   /*
   * Destructor.
   */
   ModifierManager::~ModifierManager()
   {} 

   /*
   * Read parameter file. 
   *
   * \param in input parameter file stream.
   */
   void ModifierManager::readParam(std::istream &in)
   {
      readBegin(in, "ModifierManager");
      Manager<Modifier>::readParam(in);
   }

   // Setup

   void ModifierManager::setupPostExchange() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].setupPostExchange();
      }
   }
 
   void ModifierManager::setupPostNeighbor();
   void ModifierManager::setupPostForce();
   void ModifierManager::setupEnd();

   // Integration

   void ModifierManager::sample(long iStep) 
   {
      for (int i=0; i < size(); ++i) {
         if ((*this)[i].isAtInterval(iStep)) {
            (*this)[i].preIntegrate();
         }
      } 
   }
 
   void ModifierManager::preIntegrate(long iStep);
   void ModifierManager::postIntegrate(long iStep);

   void ModifierManager::preTransform(long iStep);
   void ModifierManager::preExchange(long iStep);
   void ModifierManager::postExchange(long iStep);
   void ModifierManager::postNeighbor(long iStep);

   void ModifierManager::preUpdate(long iStep);
   void ModifierManager::postUpdate(long iStep);

   void ModifierManager::preForce(long iStep);
   void ModifierManager::postForce(long iStep);
   void ModifierManager::endOfStep(long iStep);

   void ModifierManager::postRun() 
   {
      for (int i=0; i < size(); ++i) {
         (*this)[i].postRun();
      }
   }

   // Communication
   void ModifierManager::pack_exchange(Buffer& buffer);
   void ModifierManager::unpack_exchange(Buffer& buffer);
   void ModifierManager::pack_update(Buffer& buffer);
   void ModifierManager::unpack_update(Buffer& buffer);
   void ModifierManager::pack_reverseUpdate(Buffer& buffer);
   void ModifierManager::unpack_reverseUpdate(Buffer& buffer);

   /*
   * Return pointer to default factory.
   */
   Factory<Modifier>* ModifierManager::newDefaultFactory() const
   {
      return new ModifierFactory(*simulationPtr_);
   }
 
}
#endif
