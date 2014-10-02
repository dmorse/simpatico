#ifdef  MCMD_LINK
#ifndef MCMD_LINK_H
#define MCMD_LINK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Veronica Chappa and The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/chemistry/Atom.h>

namespace McMd
{

   using namespace Util;

   //class Atom;

   /**
   * A Link represents a crosslink between two Atoms.
   *
   * All non-const mutator methods of Link are private, and are called only by
   * the friend class LinkMaster. Outside of the implementation of LinkMaster, 
   * a Link is thus a read-only object.
   *
   * \ingroup McMd_Link_Module
   */
   class Link 
   {

   public:

      /**
      * Constructor.
      */
      Link();

      /**
      * Get Atom0 connected to a Link.
      */
      const Atom& atom0() const;

      /**
      * Get Atom0 connected to a Link.
      */
      Atom& atom0();

      /**
      * Get Atom1 connected to a Link.
      */
      const Atom& atom1() const;    
 
     /**
      * Get Atom1 connected to a Link.
      */
      Atom& atom1();

      /**
      * Get the typeId for this Link.
      */
      int typeId() const;
 
      /**
      * Get a permanent integer identifier for this object.
      */
      int tag() const;
 
      /**
      * Is this Link active?
      */
      bool isActive() const;  

      /**
      * Activate or deactivate the Link.
      *
      * \param isActive true (active) or false (inactive)
      */
      void setIsActive(bool isActive);

   private:
      
      /// pointers to Atoms that belong to the Link.
      Atom*  atom0Ptr_;
      Atom*  atom1Ptr_;
   
      /// Integer index for the type of Link.
      int    typeId_;

      /// Immutable identifier.
      int    tag_;

      /// If false, the Link is deactivated.
      bool   isActive_;
   
      /**
      * Add the two atoms to the Link.
      *
      * \param atom0 atom to be added.
      * \param atom1 atom to be added.
      */
      void setAtoms(Atom &atom0, Atom &atom1);
 
      /**
      * Set the typeId for this Link.
      *
      * \param typeId 
      */
      void setTypeId(int typeId);
    
      /**
      * Set a permanent identifier for this Link.
      *
      * This function must be passed a non-negative
      * tag, and can only be called once. 
      *
      * \param tag identifier
      */
      void setTag(int tag);

      /**
      * Nullify pointers, set typeId to -1 and isActive to false.
      */
      void clear();

   //friends:

      friend class LinkMaster;

   };

   // Inline method definitions

   /*
   * Constructor.
   */
   inline Link::Link()
    : atom0Ptr_(0),
      atom1Ptr_(0),
      typeId_(-1),
      isActive_(false)
   {}

   /*
   * Add the two atoms to the Link.
   */
   inline void Link::setAtoms(Atom &atom0, Atom &atom1)
   { 
      if (&atom0 != &atom1) {
         atom0Ptr_ = &atom0; 
         atom1Ptr_ = &atom1;
      }
   }

   /*
   * Get a const reference to atom0.
   */
   inline const Atom& Link::atom0() const 
   {  return *atom0Ptr_; }

   /*
   * Get atom0.
   */
   inline Atom& Link::atom0()
   {  return *atom0Ptr_; }

   /*
   * Get a const referencen to atom1.
   */
   inline const Atom& Link::atom1() const 
   {  return *atom1Ptr_; }

   /*
   * Get atom1.
   */
   inline Atom& Link::atom1()
   {  return *atom1Ptr_; }

   /*
   * Set the type id for this Link.
   */
   inline void Link::setTypeId(int typeId)
   {  typeId_ = typeId; }

   /*
   * Get the typeId for this Link.
   */ 
   inline int Link::typeId() const
   {  return typeId_; }

   /*
   * Is this Link active?
   */
   inline bool Link::isActive() const 
   {  return isActive_; }  

   /*
   * Activate or deactivate the Link.
   */
   inline void Link::setIsActive(bool isActive)
   {  isActive_ = isActive; }

   /*
   * Set a permanent identifier for this Link.
   */
   inline void Link::setTag(int tag)
   {  tag_ = tag; }

   /*
   * Get a permanent integer identifier for this object.
   */
   inline int Link::tag() const
   {  return tag_;  }

   /*
   * Nullify pointers, set typeId to -1 and isActive to false.
   */
   inline void Link::clear()
   {
      atom0Ptr_ = 0;
      atom1Ptr_ = 0;
      typeId_   =-1;
      isActive_ = false;
   }

} 
#endif
#endif
