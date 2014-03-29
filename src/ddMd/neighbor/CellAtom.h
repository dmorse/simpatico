#ifndef DDMD_CELL_ATOM_H
#define DDMD_CELL_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace DdMd
{

   /**
   * Data for an atom in a CellList.
   */
   class CellAtom
   {

   public:

      #if 0
      struct Tag {
         Atom* ptr;
         int cellRank;
      };
      #endif

      void setPtr(Atom* atomPtr) 
      {  ptr_ = atomPtr; }

      void update() 
      {
         position_ = ptr_->position();
         id_ = ptr_->id();
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

      Mask* maskPtr() const
      {  
         return &(ptr_->mask()); 
      }

   private:

      Vector position_;
      Atom* ptr_;
      int id_;

   };

}
#endif
