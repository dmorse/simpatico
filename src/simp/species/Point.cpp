/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Point.h"

namespace Simp
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
      // nBond_, nAngle_, & nDihedral_ initialized to 0 in Species.
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

}
