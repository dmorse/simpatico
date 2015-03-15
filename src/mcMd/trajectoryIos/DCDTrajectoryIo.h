#ifndef MCMD_DCD_TRAJECTORY_IO_H
#define MCMD_DCD_TRAJECTORY_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/trajectoryIos/TrajectoryIo.h> // base class
#include <util/boundary/Boundary.h>      // base class
#include <util/containers/DArray.h>

namespace McMd
{

   using namespace Util;
   
   /**
   * TrajectoryIo for CHARMM ".dcd" data files. Currently it has only been tested
   * with Hoomd generated trajectory files.
   *
   * It is assumed that the atoms belonging to all species and molecules of a species occur in the file
   * consecutively.
   *
   * \ingroup McMd_TrajectoryIo_Module
   */
   class DCDTrajectoryIo : public TrajectoryIo
   {
   
   public:

      /**
      * Constructor. 
      */
      DCDTrajectoryIo(System& system);

      /** 
      * Destructor.   
      */
      virtual ~DCDTrajectoryIo();
 
      /**
      * Read trajectory file header and initialize simulation parameters.
      *
      * \param file input file stream
      */
      void readHeader(std::fstream &file);
 
      /**
      * Read a single frame. Frames are assumed to be read consecutively.
      *
      * \param file input file stream
      * \return true if another frame exists, false if at end
      */
      bool readFrame(std::fstream &file);

      private:

          /// The number of atoms stored in the file.
          int nAtoms_;

          /// The total number of frames in the trajectory file.
          int nFrames_;

          /// The current frame index, starting from 0.
          int frameId_;

          /// buffer for reading data (x values).
          DArray<float> xBuffer_;

          /// buffer for reading data (y values).
          DArray<float> yBuffer_;

          /// buffer for reading data (z values).
          DArray<float> zBuffer_;

   }; 

} 
#endif
