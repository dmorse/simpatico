#ifndef MCMD_TETHER_H
#define MCMD_TETHER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>    // member

namespace McMd
{

   using namespace Util;

   class Atom;

   /**
   * A Tether represents an interaction between an Atom and an anchor position.
   *
   * An active Tether has an anchor position Vector and a pointer to an 
   * associated Atom. A Tether represents an interaction that depends upon the
   * distance between the atom and the anchor position, but does not define the
   * interaction potential.
   *
   * A Tether object may also optionally be associated with another Tether, its
   * partner, as done in some slip-link algorithms. 
   *
   * All non-const mutator methods of Tether are private, and are called only by
   * the friend class TetherMaster. Outside of the implementation of TetherMaster, 
   * a Tether is thus a read-only object.
   *
   * \ingroup Tether_Module
   */
   class Tether 
   {

   public:

      /**
      * Constructor.
      */
      Tether();

      /**
      * Get anchor position Vector by reference.
      */
      const Vector& anchor() const;
   
      /**
      * Get the tethered Atom by reference.
      */
      Atom& atom() const;
   
      /**
      * Get the tether type id.
      */
      int typeId() const;
   
      /**
      * Get the partner (if any) by reference.
      */
      Tether& partner() const;

      /**
      * Is this Tether associated with Atom?
      */
      bool hasAtom() const;

      /**
      * Is Tether associated with a partner Tether?
      */
      bool hasPartner() const;

   private:

      /**
      * Position of anchor.
      */
      Vector anchor_;

      /**
      * Pointer to tethered Atom.
      */
      Atom*   atomPtr_;

      /**
      * Pointer to partner Tether.
      * 
      * Used only in slip link algorithms with pairs of associated tethers.
      * If a Tether has no partner, this pointer should be null.
      */
      Tether* partnerPtr_;

      /**
      * Tether type id.
      */
      int     typeId_;

      /**
      * Set the anchor.
      *
      * \param anchor anchor position Vector.
      */ 
      void setAnchor(const Vector& anchor);

      /**
      * Associate a Tether with an Atom.
      *
      * \param atom tethered atom
      */ 
      void setAtom(Atom& atom);

      /**
      * Nullify pointer to associated Atom.
      */ 
      void unsetAtom();

      /**
      * Set the tether type Id.
      *
      * \param typeId tether type id.
      */ 
      void setTypeId(int typeId);

      /**
      * Associate a Tether with a partner.
      *
      * \param partner partner Tether object
      */ 
      void setPartner(Tether& partner);

      /**
      * Nullify pointer to a partner Tether.
      */ 
      void unsetPartner();

   //friends:

      friend class TetherMaster;

   };

   // Inline method definitions

   /*
   * Constructor.
   */
   inline Tether::Tether()
    : anchor_(),
      atomPtr_(0),
      partnerPtr_(0),
      typeId_(0)
   {}

   /*
   * Get anchor position by reference.
   */
   inline void Tether::setAnchor(const Vector& anchor) 
   {  anchor_ = anchor; }

   /*
   * Associate this Tether with an Atom.
   */
   inline void Tether::setAtom(Atom& atom)
   {  atomPtr_ = &atom; }

   /*
   * Nullify pointer to Atom.
   */ 
   inline void Tether::unsetAtom()
   {  if (atomPtr_ != 0) atomPtr_ = 0; }

   /*
   * Set Tether type index.
   */
   inline void Tether::setTypeId(int typeId)
   {  typeId_ = typeId; }

   /*
   * Associate this Tether with a partner Tether.
   */
   inline void Tether::setPartner(Tether& partner)
   {  partnerPtr_         = &partner; }

   /*
   * Nullify pointer to Partner.
   */ 
   inline void Tether::unsetPartner()
   {  if (partnerPtr_ != 0) partnerPtr_ = 0; }

   /*
   * Get anchor position by reference.
   */
   inline const Vector& Tether::anchor() const
   {  return anchor_; }

   /*
   * Get associated Atom by reference.
   */
   inline Atom& Tether::atom() const
   {  return *atomPtr_; }

   /*
   * Get type id.
   */
   inline int Tether::typeId() const
   {  return typeId_; }

   /*
   * Get partner Tether by reference.
   */
   inline Tether& Tether::partner() const
   {  return *partnerPtr_; }

   /*
   * Return true if this Tether is associated with an Atom.
   */
   inline bool Tether::hasAtom() const
   {  return (atomPtr_ != 0); }

   /*
   * Return true if this Tether has a partner.
   */
   inline bool Tether::hasPartner() const
   {  return (partnerPtr_ != 0); }


} 
#endif
