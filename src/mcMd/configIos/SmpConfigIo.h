#ifndef MCMD_SMP_CONFIG_IO_H
#define MCMD_SMP_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/configIos/ConfigIo.h>
#include <simp/boundary/Boundary.h>
#include <util/global.h>

#include <iostream>

namespace McMd
{

   class System;
   class Atom;
   
   using namespace Util;
   using namespace Simp;

   /**
   * Common configuration file format for simpatico.
   *
   * \ingroup McMd_ConfigIo_Module
   */
   class SmpConfigIo : public ConfigIo
   {
   
   public:

      /// Constructor. 
      SmpConfigIo(System& system);
 
      /// Destructor.   
      virtual ~SmpConfigIo();
 
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

   }; 

} 
#endif
