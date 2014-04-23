#ifndef MDPP_GROUP_H
#define MDPP_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace MdPp
{

   /**
   * A group of covalently interacting atoms.
   */
   template <int N>
   class Group
   {
   public: 
   
      /**   
      * Array of integer ids of atoms in this group.
      */
      int atomIds[N];
  
      /** 
      * Integer index for the type of group.
      */
      int typeId;
  
      /** 
      * Global id for this group.
      */
      int id;

   };

} 
#endif
