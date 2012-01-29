#ifndef LATTICE_SYSTEM_H
#define LATTICE_SYSTEM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

namespace DdMd
{

   /**
   * Enumeration of the 7 possible Bravais lattice systems.
   *
   * \ingroup Boundary_Module
   */
   enum LatticeSystem {Cubic, Tetragonal, Orthorhombic, 
                       Monoclinic, Triclinic, Rhombohedral, Hexagonal};


   /**
   * istream extractor for a LatticeSystem.
   *
   * \param  in       input stream
   * \param  lattice  LatticeSystem to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, LatticeSystem& lattice);

   /**
   * ostream inserter for an LatticeSystem.
   *
   * \param  out      output stream
   * \param  lattice  LatticeSystem to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, LatticeSystem lattice);

}
#endif
