#ifndef MCMD_MD_CONFIG_IO_H
#define MCMD_MD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/configIos/McMdConfigIo.h>
#include <util/global.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   class System;
   class Atom;
   
   /**
   * ConfigIo for MD simulations (includes velocities).
   *
   * \ingroup McMd_ConfigIo_Module
   */
   class MdConfigIo : public McMdConfigIo
   {
   
   public:

      /// Constructor. 
      MdConfigIo(System& system);
 
      /// Destructor.   
      virtual ~MdConfigIo();
 
   protected:

      /// Read data for one atom.
      virtual void readAtom(std::istream& out, Atom& atom);
     
      /// Write data for one atom.
      virtual void writeAtom(std::ostream& out, const Atom& atom);
     
   }; 

} 
#endif
