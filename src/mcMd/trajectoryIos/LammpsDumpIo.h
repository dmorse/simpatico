#ifndef MCMD_LAMMPS_DUMP_IO_H
#define MCMD_LAMMPS_DUMP_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/trajectoryIos/TrajectoryIo.h> // base class
#include <util/containers/DArray.h>

namespace McMd
{

   using namespace Util;

   /**
   * TrajectoryIo for Lammps dump trajectory files.
   *
   * This class assumes that atom tags are ordered by molecule and species, 
   * with consecutive ids for atoms in the same molecule and consecutive 
   * blocks for molecules in the same species.
   *
   * \ingroup McMd_TrajectoryIo_Module
   */
   class LammpsDumpIo : public TrajectoryIo
   {
   
   public:

      /**
      * Constructor. 
      */
      LammpsDumpIo(System& system);

      /** 
      * Destructor.   
      */
      virtual ~LammpsDumpIo();
 
      /**
      * Setup before reading frames.
      *
      * \param file input file stream.
      */
      void readHeader(std::fstream& file);

      /**
      * Read a single frame. Frames are assumed to be read consecutively.
      *
      * \param file input file stream.
      * \return true if a frame is available, false if end of file
      */
      bool readFrame(std::fstream& file);

   private:

       /// Atom positions, indexed by id.
       DArray< Vector > positions_;

   }; 

} 
#endif
