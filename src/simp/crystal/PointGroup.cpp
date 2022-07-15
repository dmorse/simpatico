/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PointGroup.h"

namespace Simp {

   /**
   * Output stream inserter operator writing a PointGroup to stream
   */ 
   std::ostream& operator << (std::ostream& out, const PointGroup& g)
   {
      int i, size;
      size = g.size();
      out << "size = " << size << std::endl;
      for (i = 0; i < size; ++i) {
         out << std::endl;
         out << g[i];
      }
      return out;
   }

}
