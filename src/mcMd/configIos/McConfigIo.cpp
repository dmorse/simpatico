#ifndef MCMD_MC_CONFIG_IO_CPP
#define MCMD_MC_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * Read data for one atom and transform its position from
   * cartesian to generalized system.
   */
   void McConfigIo::transformCartToGen(std::istream &in, Atom& atom)
   {  
      Vector cartPosition, genPosition;
      in >> cartPosition;
      boundary().shift(cartPosition);
      boundary().transformCartToGen(cartPosition, genPosition);
      atom.position() = genPosition;
   }

   /* 
   * Write data for one atom.
   */
   void McConfigIo::writeAtom(std::ostream &out, const Atom& atom)
   {  out << atom.position() << std::endl; }

   /* 
   * Transform atom position from generalized to cartesian
   * system and write data for atom.
   */
   void McConfigIo::transformGenToCart(std::ostream &out, const Atom& atom)
   {  
      Vector genPosition, cartPosition;
      genPosition = atom.position();
      boundary().transformGenToCart(genPosition, cartPosition);
      out << cartPosition << std::endl;
   }
} 
#endif
