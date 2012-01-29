#ifndef MC_MD_CONFIG_IO_H
#define MC_MD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/configIos/ConfigIo.h>
#include <mcMd/boundary/Boundary.h>
#include <util/global.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   class System;
   class Atom;
   
   /**
   * Base class for default Mc and Md configIos.
   *
   * \ingroup ConfigIo_Module
   */
   class McMdConfigIo : public ConfigIo
   {
   
   public:

      /// Constructor. 
      McMdConfigIo(System& system);
 
      /// Destructor.   
      virtual ~McMdConfigIo();
 
      /**
      * Read configuration (particle positions) from file.
      *
      * \param in input file stream.
      */
      virtual void read(std::istream &in);
 
      /**
      * Write configuration (particle positions) to file.
      *
      * \param out output file stream.
      */
      virtual void write(std::ostream& out);
 
   protected:

      /// Read data for one atom.
      virtual void readAtom(std::istream& out, Atom& atom) = 0;
     
      /// Write data for one atom.
      virtual void writeAtom(std::ostream& out, const Atom& atom) = 0;
     
   }; 

} 
#endif
