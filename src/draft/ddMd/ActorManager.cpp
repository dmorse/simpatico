#ifndef DDMD_ACTOR_MANAGER_CPP
#define DDMD_ACTOR_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ActorManager.h" 
#include "ActorFactory.h" 

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ActorManager::ActorManager(Simulation& simulation)
   : Manager<Actor>(),
     simulationPtr_(&simulation)
   {}

   /*
   * Destructor.
   */
   ActorManager::~ActorManager()
   {} 

   /*
   * Read parameter file. 
   *
   * \param in input parameter file stream.
   */
   void ActorManager::readParam(std::istream &in)
   {
      readBegin(in, "ActorManager");
      Manager<Actor>::readParam(in);

      Actor* ptr;
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

   void ActorManager::setupPostExchange() 
   {
      int n = actionsSetupPostExchange_.size();
      for (int i = 0; i < n; ++i) {
         actionsSetupPostExchange_[i]->setupPostExchange();
      }
   }
 
   void ActorManager::setupPostNeighbor()
   {
      int n = actionsSetupPostNeighbor_.size();
      for (int i = 0; i < n; ++i) {
         actionsSetupPostNeighbor_[i]->setupPostNeighbor();
      }
   }

   void ActorManager::setupPostForce()
   {
      int n = actionsSetupPostForce_.size();
      for (int i = 0; i < n; ++i) {
         actionsSetupPostForce_[i]->setupPostForce();
      }
   }

   // Integration

   void ActorManager::preIntegrate(long iStep)
   {
      Actor* ptr;
      int n = actionsPreIntegrate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreIntegrate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preIntegrate();
         }
      }
   }
 
   void ActorManager::postIntegrate(long iStep)
   {
      Actor* ptr;
      int n = actionsPostIntegrate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostIntegrate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postIntegrate();
         }
      }
   }
 
   void ActorManager::preTransform(long iStep)
   {
      Actor* ptr;
      int n = actionsPreTransform_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreTransform_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preTransform();
         }
      }
   }
 
   void ActorManager::preExchange(long iStep)
   {
      Actor* ptr;
      int n = actionsPreExchange_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreExchange_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preExchange();
         }
      }
   }
 
   void ActorManager::postExchange(long iStep)
   {
      Actor* ptr;
      int n = actionsPostExchange_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostExchange_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postExchange();
         }
      }
   }
 
   void ActorManager::postNeighbor(long iStep)
   {
      Actor* ptr;
      int n = actionsPostNeighbor_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostNeighbor_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postNeighbor();
         }
      }
   }
 
   void ActorManager::preUpdate(long iStep)
   {
      Actor* ptr;
      int n = actionsPreUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preUpdate();
         }
      }
   }
 
   void ActorManager::postUpdate(long iStep)
   {
      Actor* ptr;
      int n = actionsPostUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postUpdate();
         }
      }
   }
 
   void ActorManager::preForce(long iStep)
   {
      Actor* ptr;
      int n = actionsPreForce_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPreForce_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preForce();
         }
      }
   }
 
   void ActorManager::postForce(long iStep)
   {
      Actor* ptr;
      int n = actionsPostForce_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPostForce_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postForce();
         }
      }
   }
 
   void ActorManager::endOfStep(long iStep)
   {
      Actor* ptr;
      int n = actionsEndOfStep_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsEndOfStep_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->endOfStep();
         }
      }
   }
 

   // Communication

   void ActorManager::packExchange(long iStep)
   {
      Actor* ptr;
      int n = actionsPackExchange_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPackExchange_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packExchange();
         }
      }
   }
 
   void ActorManager::unpackExchange(long iStep)
   {
      Actor* ptr;
      int n = actionsUnpackExchange_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsUnpackExchange_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->unpackExchange();
         }
      }
   }
 
   void ActorManager::packUpdate(long iStep)
   {
      Actor* ptr;
      int n = actionsPackUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPackUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packUpdate();
         }
      }
   }
 
   void ActorManager::unpackUpdate(long iStep)
   {
      Actor* ptr;
      int n = actionsUnpackUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsUnpackUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->unpackUpdate();
         }
      }
   }

   void ActorManager::packReverseUpdate(long iStep)
   {
      Actor* ptr;
      int n = actionsPackReverseUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsPackReverseUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packReverseUpdate();
         }
      }
   }
 
   void ActorManager::unpackReverseUpdate(long iStep)
   {
      Actor* ptr;
      int n = actionsUnpackReverseUpdate_.size();
      for (int i = 0; i < n; ++i) {
         ptr = actionsUnpackReverseUpdate_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->unpackReverseUpdate();
         }
      }
   }

   /*
   * Return pointer to default Actor factory.
   */
   Factory<Actor>* ActorManager::newDefaultFactory() const
   {  return new ActorFactory(*simulationPtr_); }
 
}
#endif
