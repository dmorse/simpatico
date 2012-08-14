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
      for  (int i = 0; i < size(); ++i) {
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
         if (ptr->hasPackReverseUpdate()) { 
            actionsUnpackReverseUpdate_.push_back(ptr); 
         }
      } // end for i
   }

   // Setup

   void ActionManager::setupPostExchange() 
   {
      int n = actionsSetupPostExchange_.size();
      for (int i = 0; i < n; ++i) {
         actionsSetupPostExchange_[i]->setupPostExchange();
      }
   }
 
   void ActionManager::setupPostNeighbor()
   {
      int n = actionsSetupPostNeighbor_.size();
      for (int i = 0; i < n; ++i) {
         actionsSetupPostNeighbor_[i]->setupPostNeighbor();
      }
   }

   void ActionManager::setupPostForce()
   {
      int n = actionsSetupPostForce_.size();
      for (int i = 0; i < n; ++i) {
         actionsSetupPostForce_[i]->setupPostForce();
      }
   }

   // Integration

   void ActionManager::preIntegrate(long iStep)
   {
      Action* ptr;
      int n = actionsPreIntegrate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreIntegrate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preIntegrate();
         }
      }
   }
 
   void ActionManager::postIntegrate(long iStep)
   {
      Action* ptr;
      int n = actionsPostIntegrate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostIntegrate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postIntegrate();
         }
      }
   }
 
   void ActionManager::preTransform(long iStep)
   {
      Action* ptr;
      int n = actionsPreTransform_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreTransform_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preTransform();
         }
      }
   }
 
   void ActionManager::preExchange(long iStep)
   {
      Action* ptr;
      int n = actionsPreExchange_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreExchange_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preExchange();
         }
      }
   }
 
   void ActionManager::postExchange(long iStep)
   {
      Action* ptr;
      int n = actionsPostExchange_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostExchange_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postExchange();
         }
      }
   }
 
   void ActionManager::postNeighbor(long iStep)
   {
      Action* ptr;
      int n = actionsPostNeighbor_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostNeighbor_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postNeighbor();
         }
      }
   }
 
   void ActionManager::preUpdate(long iStep)
   {
      Action* ptr;
      int n = actionsPreUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preUpdate();
         }
      }
   }
 
   void ActionManager::postUpdate(long iStep)
   {
      Action* ptr;
      int n = actionsPostUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postUpdate();
         }
      }
   }
 
   void ActionManager::preForce(long iStep)
   {
      Action* ptr;
      int n = actionsPreForce_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreForce_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preForce();
         }
      }
   }
 
   void ActionManager::postForce(long iStep)
   {
      Action* ptr;
      int n = actionsPostForce_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostForce_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postForce();
         }
      }
   }
 
   void ActionManager::endOfStep(long iStep)
   {
      Action* ptr;
      int n = actionsEndOfStep_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsEndOfStep_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->endOfStep();
         }
      }
   }
 

   // Communication

   void ActionManager::packExchange(Buffer& buffer, long iStep)
   {
      Action* ptr;
      int n = actionsPackExchange_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPackExchange_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packExchange(buffer);
         }
      }
   }
 
   void ActionManager::unpackExchange(Buffer& buffer, long iStep)
   {
      Action* ptr;
      int n = actionsUnpackExchange_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsUnpackExchange_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->unpackExchange(buffer);
         }
      }
   }
 
   void ActionManager::packUpdate(Buffer& buffer, long iStep)
   {
      Action* ptr;
      int n = actionsPackUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPackUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packUpdate(buffer);
         }
      }
   }
 
   void ActionManager::unpackUpdate(Buffer& buffer, long iStep)
   {
      Action* ptr;
      int n = actionsUnpackUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsUnpackUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->unpackUpdate(buffer);
         }
      }
   }

   void ActionManager::packReverseUpdate(Buffer& buffer, long iStep)
   {
      Action* ptr;
      int n = actionsPackReverseUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPackReverseUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packReverseUpdate(buffer);
         }
      }
   }
 
   void ActionManager::unpackReverseUpdate(Buffer& buffer, long iStep)
   {
      Action* ptr;
      int n = actionsUnpackReverseUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsUnpackReverseUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->unpackReverseUpdate(buffer);
         }
      }
   }

   #if 0
   /*
   * Return pointer to default factory.
   */
   Factory<Action>* ActionManager::newDefaultFactory() const
   {
      return new ActionFactory(*simulationPtr_);
   }
   #endif
 
}
#endif
