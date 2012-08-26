#ifndef MCMD_MC_CONFIG_IO_H
#define MCMD_MC_CONFIG_IO_H

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
   * ConfigIo for MC simulations (no atom velocities).
   *
   * \ingroup McMd_ConfigIo_Module
   */
   class McConfigIo : public McMdConfigIo
   {
   
   public:

      /// Constructor. 
      McConfigIo(System& system);
 
      /// Destructor.   
      virtual ~McConfigIo();
 
   protected:

      /// Read data for one atom.
      virtual void readAtom(std::istream& out, Atom& atom);
     
      /*
      * Read data for one atom and transform its position from 
      * cartesian to generalized system.
      */
      virtual void transformCartToGen(std::istream& out, Atom& atom);

      /// Write data for one atom.
      virtual void writeAtom(std::ostream& out, const Atom& atom);
     
      /*
      * Transform atom position from generalized to cartesian
      * system and write data for atom.
      */
      virtual void transformGenToCart(std::ostream& out, const Atom& atom);
   }; 

} 
#endif
