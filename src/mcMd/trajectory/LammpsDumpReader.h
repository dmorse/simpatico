#ifndef MCMD_LAMMPS_DUMP_READER_H
#define MCMD_LAMMPS_DUMP_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/trajectory/TrajectoryReader.h> // base class
#include <util/containers/DArray.h>           // member template
#include <util/space/Vector.h>                // template argument

namespace McMd
{

   using namespace Util;

   /**
   * TrajectoryReader for Lammps dump trajectory files.
   *
   * This class assumes that atom ids are ordered by molecule and species, 
   * with consecutive ids for atoms in the same molecule and consecutive 
   * blocks for molecules in the same species. It does not require that 
   * the atoms be ordered consecutively by id within the dump file.
   *
   * \ingroup McMd_TrajectoryReader_Module
   */
   class LammpsDumpReader : public TrajectoryReader
   {
   
   public:

      /**
      * Constructor. 
      */
      LammpsDumpReader(System& system);

      /** 
      * Destructor.   
      */
      virtual ~LammpsDumpReader();
 
      /**
      * Open file.
      *
      * \param filename input file name.
      */
      void open(std::string filename);

      /**
      * Read a single frame. Frames are assumed to be read consecutively.
      *
      * \return true if this frame is available, false if end of file
      */
      bool readFrame();

      /**
      * Close file.
      */
      void close();

   private:

       /// Trajectory file.
       std::ifstream file_;

       /// Atom positions, indexed by id.
       DArray< Vector > positions_;

   }; 

} 
#endif
