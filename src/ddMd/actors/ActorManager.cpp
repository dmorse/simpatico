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
   * Default constructor.
   */
   ActorManager::ActorManager()
   : Manager<Actor>(),
     simulationPtr_(0),
     setupPostExchangeActors_(),
     setupPostNeighborActors_(),
     setupPostForceActors_(),
     preIntegrate1Actors_(),
     postIntegrate1Actors_(),
     preTransformActors_(),
     preExchangeActors_(),
     postExchangeActors_(),
     postNeighborActors_(),
     preUpdateActors_(),
     postUpdateActors_(),
     preForceActors_(),
     postForceActors_(),
     endOfStepActors_(),
     exchangeActors_(),
     updateActors_(),
     reverseUpdateActors_()
   {}

   /*
   * Constructor.
   */
   ActorManager::ActorManager(Simulation& simulation)
   : Manager<Actor>(),
     simulationPtr_(&simulation),
     setupPostExchangeActors_(),
     setupPostNeighborActors_(),
     setupPostForceActors_(),
     preIntegrate1Actors_(),
     postIntegrate1Actors_(),
     preTransformActors_(),
     preExchangeActors_(),
     postExchangeActors_(),
     postNeighborActors_(),
     preUpdateActors_(),
     postUpdateActors_(),
     preForceActors_(),
     postForceActors_(),
     endOfStepActors_(),
     exchangeActors_(),
     updateActors_(),
     reverseUpdateActors_()
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
      Manager<Actor>::readParameters(in);

      Actor* ptr;
      for  (int i = 0; i < size(); ++i) {
         ptr = &(*this)[i];
         if (ptr->isSet(Actor::Flags::SetupPostExchange)) { 
            setupPostExchangeActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::SetupPostNeighbor)) { 
            setupPostNeighborActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::SetupPostForce)) { 
            setupPostForceActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PreIntegrate1)) { 
            preIntegrate1Actors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PostIntegrate1)) { 
            postIntegrate1Actors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PreTransform)) { 
            preTransformActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PreExchange)) { 
            preExchangeActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PostExchange)) { 
            postExchangeActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PostNeighbor)) { 
            postNeighborActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PreUpdate)) { 
            preUpdateActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PostUpdate)) { 
            postUpdateActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PreForce)) { 
            preForceActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::PostForce)) { 
            postForceActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::EndOfStep)) { 
            endOfStepActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::Exchange)) { 
            exchangeActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::Update)) { 
            updateActors_.append(*ptr); 
         }
         if (ptr->isSet(Actor::Flags::ReverseUpdate)) { 
            reverseUpdateActors_.append(*ptr); 
         }
      } // end for i
   }

   // Setup

   void ActorManager::setupPostExchange() 
   {
      int n = setupPostExchangeActors_.size();
      for (int i = 0; i < n; ++i) {
         setupPostExchangeActors_[i].setupPostExchange();
      }
   }
 
   void ActorManager::setupPostNeighbor()
   {
      int n = setupPostNeighborActors_.size();
      for (int i = 0; i < n; ++i) {
         setupPostNeighborActors_[i].setupPostNeighbor();
      }
   }

   void ActorManager::setupPostForce()
   {
      int n = setupPostForceActors_.size();
      for (int i = 0; i < n; ++i) {
         setupPostForceActors_[i].setupPostForce();
      }
   }

   // Integration

   void ActorManager::preIntegrate1(long iStep)
   {
      Actor* ptr;
      int n = preIntegrate1Actors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preIntegrate1Actors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preIntegrate1();
         }
      }
   }
 
   void ActorManager::postIntegrate1(long iStep)
   {
      Actor* ptr;
      int n = postIntegrate1Actors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postIntegrate1Actors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postIntegrate1();
         }
      }
   }
 
   void ActorManager::preTransform(long iStep)
   {
      Actor* ptr;
      int n = preTransformActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preTransformActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preTransform();
         }
      }
   }
 
   void ActorManager::preExchange(long iStep)
   {
      Actor* ptr;
      int n = preExchangeActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preExchangeActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preExchange();
         }
      }
   }
 
   void ActorManager::postExchange(long iStep)
   {
      Actor* ptr;
      int n = postExchangeActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postExchangeActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postExchange();
         }
      }
   }
 
   void ActorManager::postNeighbor(long iStep)
   {
      Actor* ptr;
      int n = postNeighborActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postNeighborActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postNeighbor();
         }
      }
   }
 
   void ActorManager::preUpdate(long iStep)
   {
      Actor* ptr;
      int n = preUpdateActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preUpdateActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preUpdate();
         }
      }
   }
 
   void ActorManager::postUpdate(long iStep)
   {
      Actor* ptr;
      int n = postUpdateActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postUpdateActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postUpdate();
         }
      }
   }
 
   void ActorManager::preForce(long iStep)
   {
      Actor* ptr;
      int n = preForceActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preForceActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preForce();
         }
      }
   }
 
   void ActorManager::postForce(long iStep)
   {
      Actor* ptr;
      int n = postForceActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postForceActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postForce();
         }
      }
   }
 
   void ActorManager::endOfStep(long iStep)
   {
      Actor* ptr;
      int n = endOfStepActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &endOfStepActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->endOfStep();
         }
      }
   }

   // Communication

   void ActorManager::packExchange(long iStep)
   {
      Actor* ptr;
      int n = exchangeActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &exchangeActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packExchange();
         }
      }
   }
 
   void ActorManager::unpackExchange(long iStep)
   {
      Actor* ptr;
      int n = exchangeActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &exchangeActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->unpackExchange();
         }
      }
   }
 
   void ActorManager::packUpdate(long iStep)
   {
      Actor* ptr;
      int n = updateActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &updateActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packUpdate();
         }
      }
   }
 
   void ActorManager::unpackUpdate(long iStep)
   {
      Actor* ptr;
      int n = updateActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &updateActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->unpackUpdate();
         }
      }
   }

   void ActorManager::packReverseUpdate(long iStep)
   {
      Actor* ptr;
      int n = reverseUpdateActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &reverseUpdateActors_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->packReverseUpdate();
         }
      }
   }
 
   void ActorManager::unpackReverseUpdate(long iStep)
   {
      Actor* ptr;
      int n = reverseUpdateActors_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &reverseUpdateActors_[i];
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
