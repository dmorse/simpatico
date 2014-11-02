#ifdef  MCMD_LINK
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Veronica Chappa and The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinkMaster.h"         
#include  <util/misc/Observer.h>
#include  <util/random/Random.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   LinkMaster::LinkMaster()
    : linkCapacity_(0),
      atomCapacity_(0)
   { setClassName("LinkMaster"); }

   /**
   * Read linkCapacity and allocate.
   */
   void LinkMaster::readParameters(std::istream& in)
   {
      read<int>(in, "linkCapacity", linkCapacity_); 
      read<int>(in, "atomCapacity", atomCapacity_); 
      allocate();
   }

   /**
   * Add a new link.
   */ 
   void LinkMaster::addLink(Atom& atom0, Atom& atom1, int typeId)
   {
      Link* linkPtr;
      int   atom0Id = atom0.id();
      int   atom1Id = atom1.id();

      // Check preconditions
      if (atom0Id < 0 || atom0Id >= atomCapacity_) {
         Log::file() << "Atom0Id       = " << atom0Id << std::endl;
         Log::file() << "atomCapacity = " << atomCapacity_ << std::endl;
         UTIL_THROW("Invalid atom0 id");
      }
      if (atom1Id < 0 || atom1Id >= atomCapacity_) {
         Log::file() << "Atom1Id       = " << atom1Id << std::endl;
         Log::file() << "atomCapacity = " << atomCapacity_ << std::endl;
         UTIL_THROW("Invalid atom1 id");
      }

      // Pop an unused Link off the reservoir
      linkPtr = &reservoir_.pop();
      linkSet_.append(*linkPtr);

      // Initialize the link, and add its address to linkPtrs_[atomId]
      linkPtr->setAtoms(atom0,atom1);
      linkPtr->setTypeId(typeId);
      linkPtr->setIsActive(true);
      atomLinkSets_[atom0Id].append(*linkPtr);
      atomLinkSets_[atom1Id].append(*linkPtr);

      // Notify observers of addition of this Link.
      LinkAddEvent event(linkPtr);
      Notifier<LinkAddEvent>::notifyObservers(event);
   }


   /**
   * Remove a Link.
   */
   void LinkMaster::removeLink(int id)
   {
      Link* linkPtr = &link(id);
      int  atom0Id;
      int  atom1Id;

      if (!link(id).isActive()) {
        UTIL_THROW("Attempt to remove a nonactive link");
      }
      atom0Id = linkPtr->atom0().id();
      atom1Id = linkPtr->atom1().id();
      // Remove link from atom0 and atom1 link sets
      if (!atomLinkSets_[atom0Id].isElement(*linkPtr)) {
        UTIL_THROW("Link is not in atomLinkSets of atom0");
      }
      if (!atomLinkSets_[atom1Id].isElement(*linkPtr)) {
        UTIL_THROW("Link is not in atomLinkSets of atom1");
      }
      atomLinkSets_[atom0Id].remove(*linkPtr);
      atomLinkSets_[atom1Id].remove(*linkPtr);

      // Notify observers of removal of this Link.
      // Notify before clearing Link so atomPtrs and typeId are available.
      LinkRemoveEvent event(linkPtr);
      Notifier<LinkRemoveEvent>::notifyObservers(event);

      // Clear the link: nullify atomPtrs, set typeId = -1, isActive = false
      linkPtr->clear();

      // Return link to reservoir 
      linkSet_.remove(*linkPtr);
      reservoir_.push(*linkPtr);

   }

   /*
   * Modify the atoms attached to a link
   */
   void LinkMaster::reSetAtoms(Link& link, Atom& atom0, Atom& atom1)
   {
      int   atom0Id = atom0.id();
      int   atom1Id = atom1.id();
      int   oldAtom0Id = link.atom0().id();
      int   oldAtom1Id = link.atom1().id();
     
      // Check preconditions    
      if (!link.isActive()) {
        UTIL_THROW("Attempt to reset atoms for a nonactive link");
      } 
      if (atom0Id < 0 || atom0Id >= atomCapacity_) {
         Log::file() << "Atom0Id       = " << atom0Id << std::endl;
         Log::file() << "atomCapacity = " << atomCapacity_ << std::endl;
         UTIL_THROW("Invalid atom0 id");
      }
      if (atom1Id < 0 || atom1Id >= atomCapacity_) {
         Log::file() << "Atom1Id       = " << atom1Id << std::endl;
         Log::file() << "atomCapacity = " << atomCapacity_ << std::endl;
         UTIL_THROW("Invalid atom1 id");
      }

      // Remove link from oldAtom0 and oldAtom1 link sets
      if (!atomLinkSets_[oldAtom0Id].isElement(link)) {
        UTIL_THROW("Link is not in atomLinkSets of atom0");
      }
      if (!atomLinkSets_[oldAtom1Id].isElement(link)) {
        UTIL_THROW("Link is not in atomLinkSets of atom1");
      }
      atomLinkSets_[oldAtom0Id].remove(link);
      atomLinkSets_[oldAtom1Id].remove(link);       

      // Change the atoms      
      link.setAtoms(atom0,atom1); 
      atomLinkSets_[atom0Id].append(link);
      atomLinkSets_[atom1Id].append(link);
      
      LinkResetEvent event(&link);
      Notifier<LinkResetEvent>::notifyObservers(event);      
 
   }
   
  /*
   * Modify one atom attached to a link
   */
   void LinkMaster::reSetAtom(Link& link, Atom& atom, int endId)
   {     
      int   atomId = atom.id();
      int   oldAtomId;
      
      if (endId==0){
        oldAtomId = link.atom0().id();
      }
      else {
        oldAtomId = link.atom1().id();
      }
      
      // Check preconditions    
      if (!link.isActive()) {
        UTIL_THROW("Attempt to reset an atom for a nonactive link");
      } 
      if (atomId < 0 || atomId >= atomCapacity_) {
         Log::file() << "AtomId       = " << atomId << std::endl;
         Log::file() << "atomCapacity = " << atomCapacity_ << std::endl;
         UTIL_THROW("Invalid atom id");
      }
 
      // Remove link from oldAtom link sets
      if (!atomLinkSets_[oldAtomId].isElement(link)) {
        UTIL_THROW("Link is not in atomLinkSets of atom");
      }     
      atomLinkSets_[oldAtomId].remove(link);       

      // Change the atom      
      if (endId==0){      
        link.setAtoms(atom,link.atom1());
      }
      else {
	link.setAtoms(link.atom0(),atom);
      }
      atomLinkSets_[atomId].append(link);
      
      ReSetAtomEvent event(&link, endId);
      Notifier<ReSetAtomEvent>::notifyObservers(event);     
 
   }      
   
   /* 
   * Get a randomly chosen link of a specified Species.
   */
   Link& LinkMaster::randomLink(Random& random)
   {
      int n = nLink(); 
      if (n <= 0) {
         UTIL_THROW("Number of links in species <= 0");
      }
      int id = random.uniformInt(0, n);
      return linkSet_[id]; 
   }


   /**
   * Allocate and initialize all required memory.
   */
   void LinkMaster::allocate() 
   {

      links_.allocate(linkCapacity_);
      reservoir_.allocate(linkCapacity_);
      linkSet_.allocate(links_);

      // Set tags for all Links.
      for (int i = 0; i < linkCapacity_; ++i) {
         links_[i].setTag(i);
      }

      // Push all links onto reservoir stack, in reverse order.
      // The first element links_[0] wll be at the top of stack.
      int i;
      for (i = linkCapacity_ - 1; i >= 0; --i) {
         reservoir_.push(links_[i]);
      }

      // Allocate array of pointers to links.
      atomLinkSets_.allocate(atomCapacity_); 
   }

   /*
   * Return true if this LinkMaster is valid, or throw an Exception.
   */
   bool LinkMaster::isValid() const 
   {

      if (links_.isAllocated()) {

         if (links_.capacity() != reservoir_.capacity()) {
            UTIL_THROW("Unequal capacities for links_ and reservoir_");
         }
         if (links_.capacity() != reservoir_.size() + linkSet_.size()) {
            UTIL_THROW("links_ capacity != reservoir_ + linkSet_ sizes");
         }

         // Check consistency of all links in linkSet
         const Link*  linkPtr;
         const Atom*  atom0Ptr;
         const Atom*  atom1Ptr;
         int   i;
         for (i = 0; i < linkSet_.size(); ++i) {

            linkPtr = &(linkSet_[i]); 
            atom0Ptr = &(linkPtr->atom0());
            atom1Ptr = &(linkPtr->atom1());

            // Check: 

            // that neither atom0Ptr nor atom1Ptr is null
            if (atom0Ptr == 0 || atom1Ptr == 0) {
               UTIL_THROW("Link with one or both atoms missing");
            }

            // that atom0Ptr and atom1Ptr are distinct
            if (atom0Ptr == atom1Ptr) {
               UTIL_THROW("Link connects an atom to itself");
            }

            // that the link is active
            if (!linkPtr->isActive()) {
               UTIL_THROW("Link in linkSet is not Active");
            }

            // that the Link is in the AtomLinkSet of atom0 and of atom1
            if (!atomLinkSets_[atom0Ptr->id()].isElement(*linkPtr)) {
               UTIL_THROW("Link is not in atomLinkSets of atom0");
            }
            if (!atomLinkSets_[atom1Ptr->id()].isElement(*linkPtr)) {
               UTIL_THROW("Link is not in atomLinkSets of atom1");
            }

         }

         // Count all active links in links_ array
         int nCount = 0;
         for (i = 0; i < linkCapacity_; ++i) {
            linkPtr = &links_[i]; 
            if (linkPtr->isActive()) {
               ++nCount;
            }
         }
         if (nCount != linkSet_.size()) {
            UTIL_THROW("Number of Links with atoms != linkSet.size()");
         }

         // Count all links int AtomLinkSets_ (each appears twice)
         nCount = 0;
         for (i = 0; i < atomCapacity_; ++i) {
            nCount += atomLinkSets_[i].size();
         }
         if (nCount != 2*linkSet_.size()) {
            UTIL_THROW("Number links in atom sets != 2*linkSet.size()");
         }

      }

      return true;
   }

   /**
   * Clear LinkMaster
   */
   void LinkMaster::clear() 
   {
      for (int id = nLink()-1; id >= 0; id--) {
         removeLink(id);
      }
   }

} 
#endif
