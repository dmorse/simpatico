#ifndef MCMD_MD_CONFIG_IO_CPP
#define MCMD_MD_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdConfigIo.h"
#include <mcMd/chemistry/Atom.h>
#include <util/space/Dimension.h>
#include <util/format/Int.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor. 
   */
   MdConfigIo::MdConfigIo(System& system)
    : McMdConfigIo(system)
   {}

   /* 
   * Destructor.   
   */
   MdConfigIo::~MdConfigIo() 
   {}
 
   /* 
   * Read data for one atom.
   */
   void MdConfigIo::readAtom(std::istream &in, Atom& atom)
   {  
      in >> atom.position();
      in >> atom.velocity();
      #ifdef MCMD_SHIFT
      in >> atom.shift();
      boundary().shift(atom.position(), atom.shift());
      #else
      boundary().shift(atom.position());
      #endif
   }

   /* 
   * Write data for one atom.
   */
   void MdConfigIo::writeAtom(std::ostream &out, const Atom& atom)
   {
      out << atom.position() << std::endl; 
      out << atom.velocity(); 
      #ifdef MCMD_SHIFT
      for (int i = 0; i < Dimension; ++i) {
         out << Int(atom.shift()[i], 3);
      }
      #endif
      out << std::endl; 
   }

} 
#endif
