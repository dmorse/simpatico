#ifndef DDMD_SP_CELL_ATOM_H
#define DDMD_SP_CELL_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/sp/chemistry/SpAtom.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace DdMd
{

   /**
   * Data for an atom in a SpCellList.
   */
   class SpCellAtom
   {

   public:

      void setPtr(SpAtom* atomPtr) 
      {  ptr_ = atomPtr; }

      void update() 
      {
         position_ = ptr_->position;
         id_ = ptr_->id;
      }

      SpAtom* ptr() const
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
      SpAtom* ptr_;
      int id_;

   };

}
#endif
