#ifndef MCMD_PMC_CONFIG_IO_H
#define MCMD_PMC_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/configIos/ConfigIo.h>  // base class
#include <util/boundary/Boundary.h>   // typedef

namespace Util {
   class Vector; 
   class OrthorhombicBoundary;
}

namespace McMd
{

   using namespace Util;

   class Simulation;
   class System;
   
   /**
   * ConfigIo for PMC configuration files.
   *
   * \ingroup McMd_ConfigIo_Module
   */
   class PmcConfigIo : public ConfigIo
   {
   
   public:

      /// Constructor. 
      PmcConfigIo(System& system);
 
      /// Destructor.   
      virtual ~PmcConfigIo();
 
      /**
      * Read configuration (particle positions) from file.
      *
      * \param in input file stream.
      */
      void read(std::istream &in);
 
      /**
      * Write configuration (particle positions) to file.
      *
      * \param out output file stream.
      */
      void write(std::ostream &out);

   }; 
 
} 
#endif
