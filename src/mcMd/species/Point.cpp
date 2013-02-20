#ifndef MCMD_POINT_CPP
#define MCMD_POINT_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * Read type.
   */
   void Point::loadSpeciesParam(Serializable::IArchive &ar)
   {  
      loadParameter<int>(ar,"type", type_); 
      allocate();
      atomTypeIds_[0] = type_;
   }
   
   /* 
   * Save atom type.
   */
   void Point::save(Serializable::OArchive &ar)
   {
      ar << moleculeCapacity_;
      ar << type_;
   }
    
   int Point::getAtomTypeId(Molecule& molecule, int index)
   {
      assert(index == 0);  
      return type_; 
   }
   
} 
#endif
