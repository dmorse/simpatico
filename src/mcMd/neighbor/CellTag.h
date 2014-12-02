#ifndef MCMD_CELL_TAG_H
#define MCMD_CELL_TAG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   /**
   * Location of the pointer to a particular Atom in a CellList.
   *
   * A CellList contains a array of CellTag objects, which are indexed by 
   * Atom Id, and an array of Cell objects, each of which stores pointers
   * of all Atoms in one cell. Each CellTag contains the information 
   * necessary to locate a particular Atom within the CellList, to allow
   * fast deletion of a specified Atom from its Cell.
   *
   * The CellTag class is used only by the CellList class.
   *
   * \ingroup McMd_Neighbor_Module
   */
   class CellTag 
   {
   
   public:

      /// Null (unknown) value for  cellId and cellPos.
      static const int NullIndex = -1;
   
      /// Cell index of Cell containing associated Atom.
      int cellId;   

      /// Array index of Atom pointer in a Cell::part_ array.
      int cellPos;  
   
      /// Default constructor.
      CellTag() 
      { clear(); }
     
      /// Set CellTag to empty/unused state.
      void clear() 
      {
         cellId  = NullIndex; 
         cellPos = NullIndex;
      }
     
   };

}
#endif
