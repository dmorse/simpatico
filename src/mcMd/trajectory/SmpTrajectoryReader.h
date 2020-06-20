#ifndef MCMD_SMP_TRAJECTORY_READER_H
#define MCMD_SMP_TRAJECTORY_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/trajectory/TrajectoryReader.h> // base class
#include <util/containers/DArray.h>           // member 
#include <util/space/Vector.h>                // member 

#include <iostream>

namespace McMd
{

   using namespace Util;

   /**
   * TrajectoryReader for a simpatico trajectory file.
   *
   * This class assumes that atom tags are ordered by molecule and species, 
   * with consecutive ids for atoms in the same molecule and consecutive 
   * blocks for molecules in the same species.
   *
   * \ingroup McMd_Trajectory_Module
   */
   class SmpTrajectoryReader : public TrajectoryReader
   {
 
   public:

      /**
      * Constructor. 
      */
      SmpTrajectoryReader(System& system);

      /** 
      * Destructor.   
      */
      virtual ~SmpTrajectoryReader();
 
      /**
      * Open trajectory file, read header, and allocate memory.
      *
      * \param filename trajectory file name
      */
      void open(std::string filename);

      /**
      * Read a single frame. Frames are assumed to be read consecutively.
      *
      * \return true if this frame is available, false if end of file
      */
      bool readFrame();

      /**
      * Close trajectory file.
      */
      void close();

   private:

      /// Trajectory file.
      std::ifstream file_;

   }; 

} 
#endif
