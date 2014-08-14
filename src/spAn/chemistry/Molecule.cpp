#ifndef SPAN_MOLECULE_CPP
#define SPAN_MOLECULE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Molecule.h"

//#define UTIL_32BIT

namespace SpAn
{

   Molecule::Molecule()
    : atoms_(0),
      speciesPtr_(0),
      id_(-1),
      nAtom_(0)
   {}
   
}
#endif
