#ifndef MCMD_DDMD_CONFIG_IO_H
#define MCMD_DDMD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/configIos/ConfigIo.h> // base class
#include <util/boundary/Boundary.h>   // base class

namespace McMd
{

   using namespace Util;

   class Simulation;
   class System;
   
   /**
   * ConfigIo for DdMd data files.
   *
   * This class reads and writes the data file format that is read
   * by the LAMMPS read_data command.  Because the DdMd data model
   * is very different than that used in Simpatico, there are 
   * signficant restrictions on the use of such files:
   *
   * 1) When reading a LAMMPS data file, the maximum total number of
   * atoms defined by the Simpatico parameter file must exactly match
   * the number of atoms in the data file.
   *
   * 2) The atom ids (or tags) in the DdMd data file must appear in 
   * order, numbered from 1.  This is checked. 
   *
   * 3) The ordering of atoms in the lammps file and the topology are
   * assumed to be consistent with that used in Simpatico. The topology
   * information in the DdMd data file is discarded. This requires
   * that atoms within each molecule be listed sequentially, with
   * molecules with the same species listed sequentially, and that
   * ordered of atoms within each molecule be consistent with the
   * numbering scheme used in Simpatico. Thus far, none of this is 
   * checked.
   *
   * \ingroup McMd_ConfigIo_Module
   */
   class DdMdConfigIo : public ConfigIo
   {
   
   public:

      /// Constructor. 
      DdMdConfigIo(System& system, bool hasMolecules = false);
 
      /// Destructor.   
      virtual ~DdMdConfigIo();
 
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

   private:

      bool hasMolecules_;

   }; 

} 
#endif
