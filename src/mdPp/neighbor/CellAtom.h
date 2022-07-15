#ifndef MDPP_CELL_ATOM_H
#define MDPP_CELL_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mdPp/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace MdPp
{

   /**
   * Data for an atom in a CellList.
   *
   * \ingroup MdPp_Neighbor_Module
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
