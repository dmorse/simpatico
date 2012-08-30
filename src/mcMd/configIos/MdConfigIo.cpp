#ifndef MCMD_MD_CONFIG_IO_CPP
#define MCMD_MD_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * Read data for one atom and transform its position from
   * cartesian to generalized system.
   */
   void MdConfigIo::transformCartToGen(std::istream &in, Atom& atom)
   {
      Vector cartPosition, genPosition;

      in >> cartPosition;
      #ifdef MCMD_SHIFT
      in >> atom.shift();
      boundary().shift(cartPosition, atom.shift());
      #else
      boundary().shift(cartPosition);
      #endif
      boundary().transformCartToGen(cartPosition, genPosition);
      atom.position() = genPosition;
      in >> atom.velocity();
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

   /* 
   * Transform atom position from generalized to cartesian
   * system and write data for atom.
   */
   void MdConfigIo::transformGenToCart(std::ostream &out, const Atom& atom)
   {
      Vector genPosition, cartPosition;
      genPosition = atom.position();
      boundary().transformGenToCart(genPosition, cartPosition);
      out << cartPosition << std::endl;
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
