#ifndef DDMD_ACTION_MANAGER_CPP
#define DDMD_ACTION_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ActionManager.h" 
//#include "ActionFactory.h" 

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ActionManager::ActionManager(Simulation& simulation)
   : Manager<Action>(),
     simulationPtr_(&simulation)
   {}

   /*
   * Destructor.
   */
   ActionManager::~ActionManager()
   {} 

   /*
   * Read parameter file. 
   *
   * \param in input parameter file stream.
   */
   void ActionManager::readParam(std::istream &in)
   {
      readBegin(in, "ActionManager");
      Manager<Action>::readParam(in);

      Action* ptr;
      if (int i = 0; i < size(); ++i) {
         ptr = &(*this)[i];
         if (ptr->hasSetupPostExchange()) { 
            actionsSetupPostExchange_.push_back(ptr); 
         }
         if (ptr->hasSetupPostNeighbor()) { 
            actionsSetupPostNeighbor_.push_back(ptr); 
         }
         if (ptr->hasSetupPostForce()) { 
            actionsSetupPostForce_.push_back(ptr); 
         }
         if (ptr->hasPreIntegrate()) { 
            actionsPreIntegrate_.push_back(ptr); 
         }
         if (ptr->hasPostIntegrate()) { 
            actionsPostIntegrate_.push_back(ptr); 
         }
         if (ptr->hasPreTransform()) { 
            actionsPreTransform_.push_back(ptr); 
         }
         if (ptr->hasPreExchange()) { 
            actionsPreExchange_.push_back(ptr); 
         }
         if (ptr->hasPostExchange()) { 
            actionsPostExchange_.push_back(ptr); 
         }
         if (ptr->hasPostNeighbor()) { 
            actionsPostNeighbor_.push_back(ptr); 
         }
         if (ptr->hasPreUpdate()) { 
            actionsPreUpdate_.push_back(ptr); 
         }
         if (ptr->hasPostUpdate()) { 
            actionsPostUpdate_.push_back(ptr); 
         }
         if (ptr->hasPreForce()) { 
            actionsPreForce_.push_back(ptr); 
         }
         if (ptr->hasPostForce()) { 
            actionsPostForce_.push_back(ptr); 
         }
         if (ptr->hasEndOfStep()) { 
            actionsEndOfStep_.push_back(ptr); 
         }
         if (ptr->hasPackExchange()) { 
            actionsUnpackExchange_.push_back(ptr); 
         }
         if (ptr->hasPackUpdate()) { 
            actionsUnpackUpdate_.push_back(ptr); 
         }
         if (ptr->hasPackreverseUpdate()) { 
            actionsUnpackreverseUpdate_.push_back(ptr); 
         }
      }
   }

   // Setup

   void ActionManager::setupPostExchange() 
   {
      Action* ptr;
      int size = actionsSetupPostExchange_.size();
      for (int i = 0; i < size; ++i) {
         actionsSetupPostExchange_[i]->setupPostExchange();
      }
   }
 
   void ActionManager::setupPostNeighbor();
   {
      Action* ptr;
      int size = actionsSetupPostNeighbor_.size();
      for (int i = 0; i < size; ++i) {
         actionsSetupPostNeighbor_[i]->setupPostNeighbor();
      }
   }

   void ActionManager::setupPostForce();
   {
      Action* ptr;
      int size = actionsSetupPostForce_.size();
      for (int i = 0; i < size; ++i) {
         actionsSetupPostForce_[i]->setupPostForce();
      }
   }

   // Integration

   void ActionManager::preIntegrate(long iStep);
   {
      Action* ptr;
      int size = actionsPreIntegrate_.size();
      for (int i = 0; i < size; ++i) {
         ptr = actionsPreIntegrate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preIntegrate();
         }
      }
   }
 
   void ActionManager::postIntegrate(long iStep){}
   void ActionManager::preTransform(long iStep){}
   void ActionManager::preExchange(long iStep){}
   void ActionManager::postExchange(long iStep){}
   void ActionManager::postNeighbor(long iStep){}
   void ActionManager::preUpdate(long iStep){}
   void ActionManager::postUpdate(long iStep){}
   void ActionManager::preForce(long iStep){}
   void ActionManager::postForce(long iStep){}
   void ActionManager::endOfStep(long iStep){}

   // Communication
   void ActionManager::packExchange(Buffer& buffer, long iStep){}
   void ActionManager::unpackExchange(Buffer& buffer, long iStep){}
   void ActionManager::packUpdate(Buffer& buffer, long iStep){}
   void ActionManager::unpackUpdate(Buffer& buffer, long iStep){}
   void ActionManager::packReverseUpdate(Buffer& buffer, long iStep){}
   void ActionManager::unpackReverseUpdate(Buffer& buffer, long iStep){}

   /*
   * Return pointer to default factory.
   */
   Factory<Action>* ActionManager::newDefaultFactory() const
   {
      return new ActionFactory(*simulationPtr_);
   }
 
}
#endif
