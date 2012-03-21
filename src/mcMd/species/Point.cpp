#ifndef POINT_CPP
#define MCMD_POINT_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Point.h"

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   Point::Point()
    : Species(),
      type_(-1)
   {  
      nAtom_ = 1; 
      nBond_ = 0;
   }
   
   /* 
   * Read type.
   */
   void Point::readSpeciesParam(std::istream &in)
   {  
      read<int>(in,"type", type_); 
      allocate();
      atomTypeIds_[0] = type_;
   }
   
   /* 
   * Return type_ for every atom.
   */
   int Point::getAtomTypeId(Molecule& molecule, int index)
   {
      assert(index == 0);  
      return type_; 
   }
   
} 
#endif
