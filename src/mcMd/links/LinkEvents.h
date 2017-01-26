#ifdef  MCMD_LINK
#ifndef MCMD_LINK_EVENTS_H
#define MCMD_LINK_EVENTS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, Veronica Chappa and The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Link;

   /**
   * Event signalling addition of Link to the LinkMaster.
   */
   class LinkAddEvent
   {

   public:

      /**
      * Constructor.
      *
      * \param ptr Pointer to new Link
      */
      LinkAddEvent(Link* ptr)
      {  ptr_ = ptr; }

      /**
      * Get pointer to newly added Link.
      */
      Link* get() const
      {  return ptr_; }

   private:

      /// Pointer to new Link.
      Link* ptr_;

   };
   
   /**
   * Event signalling reset of Link to the LinkMaster.
   */
   class LinkResetEvent
   {

   public:

      /**
      * Constructor.
      *
      * \param ptr Pointer to modified Link
      */
      LinkResetEvent(Link* ptr)
      {  ptr_ = ptr; }

      /**
      * Get pointer to modified Link.
      */
      Link* get() const
      {  return ptr_; }

   private:

      /// Pointer to modified Link.
      Link* ptr_;

   };
   
   /**
   * Event signalling removal of a Link from the LinkMaster.
   */
   struct LinkRemoveEvent 
   {

   public:

      /**
      * Constructor.
      *
      * \param ptr Pointer to deleted Link.
      */
      LinkRemoveEvent(Link* ptr)
      {  ptr_ = ptr; }

      /// Get pointer to deleted Link.
      Link* get() const
      {  return ptr_; }

   private:

      /// Pointer to deleted link.
      Link* ptr_;

   };
   
   /**
   * Event signalling the reset of an atom from the LinkMaster.
   */
   struct ReSetAtomEvent 
   {

   public:

      /**
      * Constructor.
      *
      * \param ptr   Pointer to reset link and atom.
      * \param endId end index for the reset atom.
      */
      ReSetAtomEvent(Link* ptr, int endId)
      {  ptr_ = ptr;
         Id_ = endId;
      }

      /// Get pointer to reset Link.
      Link* getLink() const
      {  return ptr_; }
      
       /// Get pointer to reset link end.
      int getEndId() const
      {  return Id_; }   

   private:

      /// Pointer to link with a reset atom.
      Link* ptr_;

      /// End id for the reset atom.
      int Id_;      

   };  

} 
#endif
#endif
