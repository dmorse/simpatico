#ifndef MCMD_MC_CONFIG_IO_CPP
#define MCMD_MC_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McConfigIo.h"
#include <mcMd/chemistry/Atom.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor. 
   */
   McConfigIo::McConfigIo(System& system)
    : McMdConfigIo(system)
   {}

   /* 
   * Destructor.   
   */
   McConfigIo::~McConfigIo() 
   {}
 
   /* 
   * Read data for one atom.
   */
   void McConfigIo::readAtom(std::istream &in, Atom& atom)
   {  
      in >> atom.position();
      boundary().shift(atom.position());
   }

   /* 
   * Write data for one atom.
   */
   void McConfigIo::writeAtom(std::ostream &out, const Atom& atom)
   {  out << atom.position() << std::endl; }

} 
#endif
