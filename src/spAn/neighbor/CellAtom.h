#ifndef SPAN_CELL_ATOM_H
#define SPAN_CELL_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace SpAn
{

   /**
   * Data for an atom in a CellList.
   *
   * \ingroup SpAn_Neighbor_Module
   */
   class CellAtom
   {

   public:

      void setPtr(Atom* atomPtr) 
      {  ptr_ = atomPtr; }

      void update() 
      {
         position_ = ptr_->position;
         id_ = ptr_->id;
      }

      Atom* ptr() const
      {  return ptr_; }

      int id() const
      {  
         //return ptr_->id(); 
         return id_;
      }

      const Vector& position() const
      {  
         //return ptr_->position(); 
         return position_;
      }

   private:

      Vector position_;
      Atom* ptr_;
      int id_;

   };

}
#endif
