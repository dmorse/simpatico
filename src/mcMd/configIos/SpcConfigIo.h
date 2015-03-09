#ifndef MCMD_SPC_CONFIG_IO_H
#define MCMD_SPC_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/configIos/ConfigIo.h>
#include <util/boundary/Boundary.h>
#include <util/global.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   class System;
   class Atom;
   
   /**
   * Base class for simpatico spc config format reader/writers.
   *
   * \ingroup McMd_ConfigIo_Module
   */
   class SpcConfigIo : public ConfigIo
   {
   
   public:

      /**
      * Constructor. 
      *
      * \param system parent System 
      */
      SpcConfigIo(System& system);

      /** 
      * Destructor.   
      */
      virtual ~SpcConfigIo();
 
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

      template <int N>
      void writeGroups(std::ostream& out);

   };

} 
#endif
