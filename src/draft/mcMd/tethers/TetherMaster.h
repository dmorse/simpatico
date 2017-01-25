#ifndef MCMD_TETHER_MASTER_H
#define MCMD_TETHER_MASTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h> 
#include "Tether.h" 
#include <util/containers/DArray.h> 
#include <util/containers/ArrayStack.h> 
#include <util/containers/ArraySet.h> 
#include <mcMd/chemistry/Atom.h> 
#include <util/global.h> 

#include <util/space/Vector.h> 

namespace McMd
{

   using namespace Util;

   class Atom;
   //class Vector;

   /**
   * TetherMaster keeps track of all Tether objects in a System.
   *
   * A TetherMaster has a private Array of Tether objects available
   * for use in a System, a provides methods to access and manage
   * these objects. It provides methods to add and remove Tethers
   * (i.e., to associate Tethers with Atoms in a system) and to
   * access active Tethers either sequentially or via the Atoms to
   * which they are attached.
   *
   * \ingroup McMd_Tether_Module
   */
   class TetherMaster : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      TetherMaster();

      // Use default destructor.

      /**
      * Read tetherCapacity and allocate.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Add a tether to a specific Atom.
      *
      * The associated Tether object is taken from a reservoir of 
      * unused tethers. 
      *
      * Preconditions: The atom must not already be attached to a tether.
      *
      * \param atom   Atom object to which Tether should be attached
      * \param anchor anchor point position Vector
      */ 
      void addTether(Atom& atom, const Vector& anchor);

      /**
      * Remove the Tether associated with an Atom.
      *
      * If the Tether has a partner, also delete the partner.
      *
      * \param atom Atom to which Tether is attached.
      */ 
      void removeTether(const Atom& atom);

      /**
      * Remove a Tether.
      *
      * If the Tether has a partner, also delete the partner.
      *
      * Precondition: The tether must be active, i.e., attached to an Atom.
      *
      * \param tether Tether object to be removed / deactivated.
      */ 
      void removeTether(Tether& tether);

      /**
      * Associate the Tethers associated with two atoms as partners.
      *
      * Precondition: Both atoms must already be tethered.
      *
      * \param atom1  first Atom object 
      * \param atom2  second Atom object 
      */ 
      void pairTethers(const Atom& atom1, const Atom& atom2);

      /**
      * Associate two Tethers as partners.
      *
      * Precondition: Both tethers must be attached to atoms.
      *
      * \param tether1  first Tether object 
      * \param tether2  second Tether object 
      */ 
      void pairTethers(Tether& tether1, Tether& tether2);

      /**
      * Transfer a tether from one Atom to another.
      *
      * Precondition:  oldAtom must be tethered, and newAtom must not.
      * Postcondition: oldAtom will be tethered, and newAtom will not.
      *
      * \param oldAtom tethered Atom on entry, untethered on return.
      * \param newAtom untethered Atom on entry, tethered on return.
      */ 
      void transferTether(Atom& oldAtom, Atom& newAtom);

      /**
      * Transfer a tether from one Atom to another.
      *
      * Precondition: tether must be attached to an atom, and newAtom must
      * be untethered.
      *
      * \param tether  Tethered object, attached to a different Atom.
      * \param newAtom untethered Atom on entry, tethered on return.
      */ 
      void transferTether(Tether& tether, Atom& newAtom);

      /**
      * Is this atom Tethered ?
      *
      * \param atom Atom object of interest.
      */ 
      bool isTethered(const Atom& atom) const;

      /**
      * Return the Tether associated with an Atom, if any.
      *
      * Precondition: The Atom must be tethered, or an Exception is thrown.
      *
      * \param atom Atom object of interest.
      */ 
      Tether& tether(const Atom& atom) const;

      /**
      * Return an active tether by an internal set index.
      *
      * An "active" Tether object is one that is attached to an Atom. 
      * Pointers to all of the active Tether objects are maintained 
      * by an ArraySet < Tether > container. At any instant, this 
      * container assigns every active tether an index in the range 
      * 0 <= id < nTether, where nTether is the total number of active 
      * tethers. The index associated with any Tether in the set can
      * change, however, whenever another Tether is removed from the
      * set by calling removeTether().
      * 
      * \param id index in the range 0 <= id < nTether.
      */ 
      Tether& tether(int id) const;

      /**
      * Get the total number of active Tethers.
      */ 
      int nTether() const;

      /**
      * Return true if this TetherMaster is valid, or throw an Exception.
      */ 
      bool isValid() const;

   private:

      /**
      * Dynamic array of Tether objects.
      */
      DArray<Tether>     tethers_;

      /**
      * Set of pointers to active Tethers in the tethers_ array.
      */
      ArraySet<Tether>   tetherSet_;

      /**
      * Stack of pointers to inactive Tethers in the tethers_ array.
      */
      ArrayStack<Tether> reservoir_;

      /**
      * Array of Tether* pointers for specific atoms, indexed by atom id.
      *
      * Element tetherPtrs_[i] points to the Tether object for atom with 
      * i = atom.id(), or is equal to 0 if the atom is untethered.
      */
      DArray<Tether*>    tetherPtrs_;

      /**
      * Allocated dimension of tethers_ and reservoir_ containers.
      */
      int tetherCapacity_;

      /**
      * Allocated dimension of tetherPtrs_ array.
      *
      * This must be equal to the dimension Atom::capacity() of the array
      * of all Atoms in the simulation. 
      */
      int atomCapacity_;

      /**
      * Allocate all required memory 
      */
      void allocate();

      /**
      * Private method to remove Tether. Called by public methods.
      */
      void removeTether(Tether* tetherPtr, int id);

   };

   // Inline methods

   /*
   * Is this atom Tethered ?
   */ 
   inline bool TetherMaster::isTethered(const Atom& atom) const
   {  return (tetherPtrs_[atom.id()]); }

   /**
   * Get the Tether associated with an Atom.
   */ 
   inline Tether& TetherMaster::tether(const Atom& atom) const
   {  
      assert(tetherPtrs_[atom.id()]);
      return (*tetherPtrs_[atom.id()]); 
   }

   /*
   * Get an active tether by index.
   */ 
   inline Tether& TetherMaster::tether(int id) const
   {  return tetherSet_[id]; }

   /*
   * Get number of active Tethers.
   */ 
   inline int TetherMaster::nTether() const
   {  return tetherSet_.size(); }

} 
#endif
