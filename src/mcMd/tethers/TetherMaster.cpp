#ifndef TETHER_MASTER_CPP
#define TETHER_MASTER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TetherMaster.h"         

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   TetherMaster::TetherMaster()
    : tetherCapacity_(0),
      atomCapacity_(0)
   {}

   /**
   * Read tetherCapacity and allocate.
   */
   void TetherMaster::readParam(std::istream& in)
   {
      readBegin(in, "TetherMaster");
      read<int>(in, "tetherCapacity", tetherCapacity_); 
      read<int>(in, "atomCapacity", atomCapacity_); 
      readEnd(in);

      allocate();
   }

   /**
   * Add a new tether.
   */ 
   void TetherMaster::addTether(Atom& atom, const Vector& anchor)
   {
      Tether* tetherPtr;
      int     atomId = atom.id();

      // Check preconditions
      if (atomId < 0 || atomId >= atomCapacity_) {
         Log::file() << "AtomId       = " << atomId << std::endl;
         Log::file() << "atomCapacity = " << atomCapacity_ << std::endl;
         UTIL_THROW("Invalid atom id");
      }
      if (tetherPtrs_[atomId]) {
         UTIL_THROW("Attempt to addTether to tethered Atom");
      }

      // Pop an unused Tether off the reservoir
      tetherPtr = &reservoir_.pop();
      tetherSet_.append(*tetherPtr);

      // Initialize the Tether, and add its address to tetherPtrs_[atomId]
      tetherPtr->setAnchor(anchor);
      tetherPtr->setAtom(atom);
      tetherPtrs_[atomId] = tetherPtr;

   }

   /**
   * Transfer a Tether from one Atom to another.
   */ 
   void TetherMaster::transferTether(Atom& oldAtom, Atom& newAtom)
   {
      Tether* tetherPtr;
      int     oldId, newId;

      // Check preconditions
      oldId = oldAtom.id(); 
      newId = newAtom.id(); 
      if (oldId < 0 || oldId >= atomCapacity_) {
         UTIL_THROW("Invalid old atom id");
      }
      if (newId < 0 || newId >= atomCapacity_) {
         UTIL_THROW("Invalid new atom id");
      }

      tetherPtr = tetherPtrs_[oldId];
      if (tetherPtr == 0) {
         UTIL_THROW("Attempt to transfer from untethered Atom");
      }
      if (tetherPtrs_[newId]) {
         UTIL_THROW("Attempt to transfer to a tethered Atom");
      }

      // Make switch
      tetherPtr->setAtom(newAtom);
      tetherPtrs_[oldId] = 0;
      tetherPtrs_[newId] = tetherPtr;

   }

   /*
   * Transfer a tether from one Atom to another.
   */ 
   void TetherMaster::transferTether(Tether& tether, Atom& newAtom)
   {
      int     oldId, newId;

      // Check preconditions
      if (!tether.hasAtom()) {
         UTIL_THROW("Attempt to transfer inactive tether");
      }
      oldId = tether.atom().id(); 
      if (oldId < 0 || oldId >= atomCapacity_) {
         UTIL_THROW("Invalid old atom id");
      }

      newId = newAtom.id(); 
      if (newId < 0 || newId >= atomCapacity_) {
         UTIL_THROW("Invalid new atom id");
      }
      if (tetherPtrs_[newId]) {
         UTIL_THROW("Attempt to transfer to a tethered Atom");
      }

      // Make switch
      tether.setAtom(newAtom);
      tetherPtrs_[oldId] = 0;
      tetherPtrs_[newId] = &tether;

   }

   /*
   * Associate two Tethered atoms as partners.
   */ 
   void TetherMaster::pairTethers(const Atom& atom1, const Atom& atom2)
   {

      // Retrieve and validate atom Ids
      int id1 = atom1.id(); 
      int id2 = atom2.id(); 
      if (id1 < 0 || id1 >= atomCapacity_) {
         UTIL_THROW("Invalid atom id for atom1");
      }
      if (id1 < 0 || id2 >= atomCapacity_) {
         UTIL_THROW("Invalid atom id for atom2");
      }

      // Retrieve and validate pointers to Tethers
      Tether* tetherPtr1 = tetherPtrs_[id1];
      Tether* tetherPtr2 = tetherPtrs_[id2];
      if (tetherPtr1 == 0) {
         UTIL_THROW("Attempt to pair untethered Atom atom1");
      }
      if (tetherPtr2 == 0) {
         UTIL_THROW("Attempt to pair untethered Atom atom2");
      }

      // Create association
      tetherPtr1->setPartner(*tetherPtr2);
      tetherPtr2->setPartner(*tetherPtr1);

   }

   /*
   * Associate two Tethers as partners.
   */ 
   void TetherMaster::pairTethers(Tether& tether1, Tether& tether2)
   {

      // Check preconditions
      if (!tether1.hasAtom()) {
         UTIL_THROW("Attempt to pair inactive tether");
      }
      if (!tether2.hasAtom()) {
         UTIL_THROW("Attempt to pair inactive tether");
      }
      if (tether1.hasPartner()) {
         UTIL_THROW("Attempt to pair Tether that has a partner");
      }
      if (tether2.hasPartner()) {
         UTIL_THROW("Attempt to pair Tether that has a partner");
      }

      // Create association
      tether1.setPartner(tether2);
      tether2.setPartner(tether1);

   }


   /**
   * Remove the Tether attached to an Atom, and its partner, if any.
   */ 
   void TetherMaster::removeTether(const Atom& atom)
   {

      Tether* tetherPtr;
      int     id;

      // Check preconditions
      id = atom.id();
      if (id < 0 || id >= atomCapacity_) {
         UTIL_THROW("Invalid atom id");
      }
      tetherPtr = tetherPtrs_[id];
      if (tetherPtr == 0) {
         UTIL_THROW("Attempt to destroy tether to untethered Atom");
      }

      // Call private remove method
      removeTether(tetherPtr, id);

   }

   /**
   * Remove a Tether, and its partner, if any.
   */ 
   void TetherMaster::removeTether(Tether& tether)
   {

      // Check preconditions
      if (!tether.hasAtom()) {
         UTIL_THROW("Attempt to destroy inactive tether");
      }
      int id = tether.atom().id();
      if (id < 0 || id >= atomCapacity_) {
         UTIL_THROW("Invalid atom id");
      }

      // Call private remove method
      removeTether(&tether, id);

   }

   /**
   * Remove a Tether, and its partner, if any.
   */ 
   void TetherMaster::removeTether(Tether* tetherPtr, int id)
   {

      // Deactivate this tether
      tetherPtrs_[id] = 0;
      tetherPtr->unsetAtom();

      // Deactivate partner, if any.
      if (tetherPtr->hasPartner()) {

         Tether* partnerPtr = &tetherPtr->partner();

	 // Break association
	 tetherPtr->unsetPartner();
         partnerPtr->unsetPartner();

	 // Deactivate partner
         id = partnerPtr->atom().id(); 
         if (id < 0 || id >= atomCapacity_) {
            UTIL_THROW("Invalid atom id");
         }
         partnerPtr->unsetAtom();
         tetherPtrs_[id] = 0;
     
	 // Return partner to reservoir 
	 tetherSet_.remove(*partnerPtr);
	 reservoir_.push(*partnerPtr);

      }

      // Return tether to reservoir 
      tetherSet_.remove(*tetherPtr);
      reservoir_.push(*tetherPtr);
   }

   /*
   * Allocate and initialize all required memory.
   */
   void TetherMaster::allocate() 
   {

      tethers_.allocate(tetherCapacity_);
      reservoir_.allocate(tetherCapacity_);
      tetherSet_.allocate(tethers_);

      // Push all tethers onto reservoir stack, in reverse order.
      // The first element tethers_[0] wll be at the top of stack.
      int i;
      for (i = tetherCapacity_ - 1; i >= 0; --i) {
         reservoir_.push(tethers_[i]);
      }

      // Allocate and nullify array of pointers to tethers.
      tetherPtrs_.allocate(atomCapacity_);
      for (i = 0; i < atomCapacity_; ++i) {
         tetherPtrs_[i] = 0;
      }

   }

   /*
   * Return true if this TetherMaster is valid, or throw an Exception.
   */
   bool TetherMaster::isValid() const
   {

      const Tether* tetherPtr;
      const Tether* partnerPtr;
      const Atom*   atomPtr;

      if (tethers_.isAllocated()) {

         if (tethers_.capacity() != reservoir_.capacity()) {
            UTIL_THROW("Unequal capacities for tethers_ and reservoir_");
         }
         if (tethers_.capacity() != reservoir_.size() + tetherSet_.size()) {
            UTIL_THROW("tethers_ capacity != reservoir_ + tetherSet_ sizes");
         }

         // Check consistency of all tethers in tetherSet
         int i, atomId;
         for (i = 0; i < tetherSet_.size(); ++i) {
            tetherPtr = &(tetherSet_[i]); 

            // Check pointer to an Atom
            if (!tetherPtr->hasAtom()) {
               UTIL_THROW("Tether in tetherSet has no Atom");
            }
            atomPtr = &tetherPtr->atom();
            atomId  = atomPtr->id();
            if (tetherPtrs_[atomId] != tetherPtr) {
               UTIL_THROW("Inconsistency a Tether and tetherPtrs_ array");
            }

            // Check partner, if any
            if (tetherPtr->hasPartner()) {
               partnerPtr = &tetherPtr->partner();
               if (!partnerPtr->hasPartner()) {
                  UTIL_THROW("Partner of Tether has no Partner");
               }
               if (tetherPtr != &partnerPtr->partner()) {
                  UTIL_THROW("Inconsistency partner addresses");
               }
            }

         }

         // Count all active tethers in tethers_ array
         int nCount = 0;
         for (i = 0; i < tetherCapacity_; ++i) {
            tetherPtr = &tethers_[i]; 
            if (tetherPtr->hasAtom()) {
               ++nCount;
            } else {
               //if (tetherPtr->hasPartner()) {
               //   UTIL_THROW("Tether has no Atom but has Partner");
               //}
            }
         }
         if (nCount != tetherSet_.size()) {
            UTIL_THROW("Number of Tethers with atoms != tetherSet.size()");
         }

         // Count all non-null pointers in tetherPtrs_ array
         nCount = 0;
         for (i = 0; i < atomCapacity_; ++i) {
            if (tetherPtrs_[i]) {
               ++nCount;
            }
         }
         if (nCount != tetherSet_.size()) {
            UTIL_THROW("Count non-null Tether* pointers != tetherSet.size()");
         }

      }
      return true;

   }

} 
#endif
