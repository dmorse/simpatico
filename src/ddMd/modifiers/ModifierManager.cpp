#ifndef DDMD_MODIFIER_MANAGER_CPP
#define DDMD_MODIFIER_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ModifierManager.h" 
#include "ModifierFactory.h" 

namespace DdMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   ModifierManager::ModifierManager()
    : Manager<Modifier>(),
      simulationPtr_(0),
      setupModifiers_(),
      preIntegrate1Modifiers_(),
      postIntegrate1Modifiers_(),
      preTransformModifiers_(),
      preExchangeModifiers_(),
      postExchangeModifiers_(),
      postNeighborModifiers_(),
      preUpdateModifiers_(),
      postUpdateModifiers_(),
      preForceModifiers_(),
      postForceModifiers_(),
      endOfStepModifiers_(),
      exchangeModifiers_(),
      updateModifiers_(),
      reverseUpdateModifiers_()
   { setClassName("ModifierManager"); }

   /*
   * Constructor.
   */
   ModifierManager::ModifierManager(Simulation& simulation)
    : Manager<Modifier>(),
      simulationPtr_(&simulation),
      setupModifiers_(),
      preIntegrate1Modifiers_(),
      postIntegrate1Modifiers_(),
      preTransformModifiers_(),
      preExchangeModifiers_(),
      postExchangeModifiers_(),
      postNeighborModifiers_(),
      preUpdateModifiers_(),
      postUpdateModifiers_(),
      preForceModifiers_(),
      postForceModifiers_(),
      endOfStepModifiers_(),
      exchangeModifiers_(),
      updateModifiers_(),
      reverseUpdateModifiers_()
   {  setClassName("ModifierManager"); }

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
   void ModifierManager::readParameters(std::istream &in)
   {
      Manager<Modifier>::readParameters(in);

      Modifier* ptr;
      for  (int i = 0; i < size(); ++i) {
         ptr = &(*this)[i];
         if (ptr->isSet(Modifier::Flags::Setup)) { 
            setupModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PreIntegrate1)) { 
            preIntegrate1Modifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PostIntegrate1)) { 
            postIntegrate1Modifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PreTransform)) { 
            preTransformModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PreExchange)) { 
            preExchangeModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PostExchange)) { 
            postExchangeModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PostNeighbor)) { 
            postNeighborModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PreUpdate)) { 
            preUpdateModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PostUpdate)) { 
            postUpdateModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PreForce)) { 
            preForceModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::PostForce)) { 
            postForceModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::EndOfStep)) { 
            endOfStepModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::Exchange)) { 
            exchangeModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::Update)) { 
            updateModifiers_.append(*ptr); 
         }
         if (ptr->isSet(Modifier::Flags::ReverseUpdate)) { 
            reverseUpdateModifiers_.append(*ptr); 
         }
      } // end for i 
   }

   // Setup

   void ModifierManager::setup() 
   {
      int n = setupModifiers_.size();
      for (int i = 0; i < n; ++i) {
         setupModifiers_[i].setup();
      }
   }
 
   // Integration actions

   void ModifierManager::preIntegrate1(long iStep)
   {
      Modifier* ptr;
      int n = preIntegrate1Modifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preIntegrate1Modifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preIntegrate1(iStep);
         }
      }
   }
 
   void ModifierManager::postIntegrate1(long iStep)
   {
      Modifier* ptr;
      int n = postIntegrate1Modifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postIntegrate1Modifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postIntegrate1(iStep);
         } 
      }
   }
 
   void ModifierManager::preTransform(long iStep)
   {
      Modifier* ptr;
      int n = preTransformModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preTransformModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preTransform(iStep);
         }
      }
   }
 
   void ModifierManager::preExchange(long iStep)
   {
      Modifier* ptr;
      int n = preExchangeModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preExchangeModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preExchange(iStep);
         }
      }
   }
 
   void ModifierManager::postExchange(long iStep)
   {
      Modifier* ptr;
      int n = postExchangeModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postExchangeModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postExchange(iStep);
         }
      }
   }
 
   void ModifierManager::postNeighbor(long iStep)
   {
      Modifier* ptr;
      int n = postNeighborModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postNeighborModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postNeighbor(iStep);
         }
      }
   }
 
   void ModifierManager::preUpdate(long iStep)
   {
      Modifier* ptr;
      int n = preUpdateModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preUpdateModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preUpdate(iStep);
         }
      }
   }
 
   void ModifierManager::postUpdate(long iStep)
   {
      Modifier* ptr;
      int n = postUpdateModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postUpdateModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postUpdate(iStep);
         }
      }
   }
 
   void ModifierManager::preForce(long iStep)
   {
      Modifier* ptr;
      int n = preForceModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &preForceModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->preForce(iStep);
         }
      }
   }
 
   void ModifierManager::postForce(long iStep)
   {
      Modifier* ptr;
      int n = postForceModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &postForceModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->postForce(iStep);
         }
      }
   }
 
   void ModifierManager::endOfStep(long iStep)
   {
      Modifier* ptr;
      int n = endOfStepModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &endOfStepModifiers_[i];
         if (ptr->isAtInterval(iStep)) {
            ptr->endOfStep(iStep);
         }
      }
   }

   // Communication

   void ModifierManager::packExchange(long iStep)
   {
      Modifier* ptr;
      int n = exchangeModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &exchangeModifiers_[i];
         ptr->packExchange();
      }
   }
 
   void ModifierManager::unpackExchange(long iStep)
   {
      Modifier* ptr;
      int n = exchangeModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &exchangeModifiers_[i];
         ptr->unpackExchange();
      }
   }
 
   void ModifierManager::packUpdate(long iStep)
   {
      Modifier* ptr;
      int n = updateModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &updateModifiers_[i];
         ptr->packUpdate();
      }
   }
 
   void ModifierManager::unpackUpdate(long iStep)
   {
      Modifier* ptr;
      int n = updateModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &updateModifiers_[i];
         ptr->unpackUpdate();
      }
   }

   void ModifierManager::packReverseUpdate(long iStep)
   {
      Modifier* ptr;
      int n = reverseUpdateModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &reverseUpdateModifiers_[i];
         ptr->packReverseUpdate();
      }
   }
 
   void ModifierManager::unpackReverseUpdate(long iStep)
   {
      Modifier* ptr;
      int n = reverseUpdateModifiers_.size();
      for (int i = 0; i < n; ++i) {
         ptr = &reverseUpdateModifiers_[i];
         ptr->unpackReverseUpdate();
      }
   }

   /*
   * Return pointer to default Modifier factory.
   */
   Factory<Modifier>* ModifierManager::newDefaultFactory() const
   {  return new ModifierFactory(*simulationPtr_); }
 
}
#endif
