/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Molecule.h"

//#define UTIL_32BIT

namespace MdPp
{

   Molecule::Molecule()
    : atoms_(0),
      speciesPtr_(0),
      id_(-1),
      nAtom_(0)
   {}
   
}
