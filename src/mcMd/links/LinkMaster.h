#ifdef  MCMD_LINK
#ifndef MCMD_LINK_MASTER_H
#define MCMD_LINK_MASTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, Veronica Chappa and David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Link.h" 
#include "LinkEvents.h" 
#include <mcMd/chemistry/Atom.h> 
#include <util/param/ParamComposite.h> 
#include <util/containers/DArray.h> 
#include <util/containers/ArrayStack.h> 
#include <util/containers/ArraySet.h> 
#include <util/containers/SSet.h> 
#include <util/containers/PArrayIterator.h> 
#include <util/containers/ConstPArrayIterator.h> 
#include <util/util/Notifier.h> 


namespace Util
{  class Random; }

namespace McMd
{

   using namespace Util;

   /**
   * Manages all Link objects in a System.
   *
   * A LinkMaster allows links between Atoms, or Link objects, to be added 
   * to or removed from a system, and provides access to existing links.
   */
   class LinkMaster : public ParamComposite, 
                      public Notifier<LinkAddEvent>,
                      public Notifier<LinkResetEvent>,                      
                      public Notifier<LinkRemoveEvent>,
                      public Notifier<ReSetAtomEvent>                     
   {

   public:

      /// Iterator for set of active links.
      typedef PArrayIterator<Link>      LinkIterator;

      /// Const iterator for set of active links.
      typedef ConstPArrayIterator<Link> ConstLinkIterator;

      /**
      * A set of links involving a particular atom.
      */
      typedef SSet<Link,400> AtomLinkSet;

      /**
      * Constructor.
      */
      LinkMaster();

      // Use default destructor.

      /**
      * Read linkCapacity and allocate.
      *
      * \param in input parameter stream
      */
      void readParam(std::istream& in);

      /**
      * Add a link betwen two specific Atoms.
      *
      * The associated Link object is taken from a reservoir of 
      * unused links. 
      *
      *
      * \param atom0  first Atom object to which a Link should be attached
      * \param atom1  second Atom object to which a Link should be attached
      * \param typeId id for the type of Link
      */ 
      void addLink(Atom& atom0, Atom& atom1, int typeId);

      /**
      * Remove a Link.
      *
      * Precondition: The Link must be active.
      *
      * \param id Link index in linkSet_, the ArraySet < Link > container.
      */ 
      void removeLink(int id);


      /**
      * Return a set of links associated with an Atom by const reference.
      *
      * \param atom Atom object of interest.
      */ 
      const AtomLinkSet& atomLinkSet(const Atom& atom) const;

      /**
      * Return a set of links associated with an Atom.
      *
      * \param atom Atom object of interest.
      */ 
      AtomLinkSet& atomLinkSet(const Atom& atom);

      /**
      * Return an active link by an internal set index.
      *
      * An "active" Link object is one that is attached to two atoms. 
      * Pointers to all of the active Link objects are maintained 
      * by an ArraySet < Link > container. At any instant, this 
      * container assigns every active link an index in the range 
      * 0 <= id < nLink, where nLink is the total number of active 
      * links. The index associated with any Link in the set can
      * change, however, whenever another Link is removed from the
      * set by calling removeLink() (not yet implemented).
      * 
      * \param id index in the range 0 <= id < nLink.
      */ 
      Link& link(int id) const;
          
      /**
      * Modify the atoms attached to a link
      */
      void reSetAtoms(Link& link, Atom& atom0, Atom& atom1);
      
      /**
      * Modify one atom attached to a link
      */
      void reSetAtom(Link& link, Atom& atom, int endId);            
            
      /**
      * Get the total number of active Links.
      */ 
      int nLink() const;

      /**
      * Initialize iterator to beginning of list of active links.
      *
      * \param iterator Iterator for active links.
      */
      void begin(LinkIterator& iterator);

      /**
      * Initialize iterator to beginning of list of active links.
      *
      * \param iterator const Iterator for active links.
      */
      void begin(ConstLinkIterator& iterator) const;

      /**
      * Get a randomly chosen link by reference.
      *
      * \param random Random number generator.
      */
      Link& randomLink(Random& random);

      /**
      * Return true if this LinkMaster is valid, or throw an Exception.
      */ 
      bool isValid() const;

      /**
      * Clear LinkMaster.
      *
      * Remove all links: clear and return them to the reservoir, empty 
      * linkSet and atomLinkSets arrays.
      */
      void clear();

   private:

      /**
      * Dynamic array of Link objects.
      */
      DArray<Link>     links_;

      /**
      * Set of pointers to active Links in the links_ array.
      */
      ArraySet<Link>   linkSet_;

      /**
      * Stack of pointers to inactive Links in the links_ array.
      */
      ArrayStack<Link> reservoir_;

      /**
      * Array of Link* pointers for specific atoms, indexed by atom id. 
      *
      * Element atomLinkSets_[i] points to the set of Links for atom with 
      * i = atom.id(), or is equal to 0 if the atom has no links.
      */
      DArray<AtomLinkSet>      atomLinkSets_;

      /**
      * Allocated dimension of links_, reservoir_ and linkSet_ containers.
      */
      int linkCapacity_;

      /**
      * Allocated dimension of atomLinkSets_ array.
      *
      * This must be equal to the dimension Atom::capacity() of the array
      * of all Atoms in the simulation. 
      */
      int atomCapacity_;

      /**
      * Allocate all required memory 
      */
      void allocate();

   };

   // Inline methods

   /*
   * Get the set of Links associated with an Atom by constant reference.
   */ 
   inline 
   const LinkMaster::AtomLinkSet& LinkMaster::atomLinkSet(const Atom& atom) const  
   {  
      return (atomLinkSets_[atom.id()]); 
   }

   /*
   * Get the set of Links associated with an Atom.
   */ 
   inline 
   LinkMaster::AtomLinkSet& LinkMaster::atomLinkSet(const Atom& atom)
   {  
      return (atomLinkSets_[atom.id()]); 
   }

   /*
   * Initialize iterator to beginning of list of active links.
   */
   inline void LinkMaster::begin(LinkMaster::LinkIterator& iterator)
   {   linkSet_.begin(iterator); }

   /**
   * Initialize iterator to beginning of list of active links.
   */
   inline void LinkMaster::begin(LinkMaster::ConstLinkIterator& iterator) const
   {   linkSet_.begin(iterator); }

   /*
   * Get an active Link by index.
   */ 
   inline Link& LinkMaster::link(int id) const
   {  return linkSet_[id]; }

   /*
   * Get number of active Links.
   */ 
   inline int LinkMaster::nLink() const
   {  return linkSet_.size(); }

} 
#endif
#endif
