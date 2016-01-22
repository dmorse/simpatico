#ifndef MCMD_DCD_TRAJECTORY_READER_H
#define MCMD_DCD_TRAJECTORY_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/trajectory/TrajectoryReader.h> // base class
#include <util/boundary/Boundary.h>           // typedef
#include <util/containers/DArray.h>           // member

#include <fstream>

namespace McMd
{

   using namespace Util;
   
   /**
   * TrajectoryReader for CHARMM ".dcd" data files. 
   *
   * It is assumed that the atoms belonging to all species and molecules of 
   * a species occur in the file consecutively. Currently it has only been 
   * tested with Hoomd generated trajectory files.
   *
   * \ingroup McMd_Trajectory_Module
   */
   class DCDTrajectoryReader : public TrajectoryReader
   {
   
   public:

      /**
      * Constructor. 
      */
      DCDTrajectoryReader(System& system);

      /** 
      * Destructor.   
      */
      virtual ~DCDTrajectoryReader();
 
      /**
      * Read trajectory file header and initialize simulation parameters.
      *
      * \param filename name of trajectory file 
      */
      void open(std::string filename);
 
      /**
      * Read a single frame. Frames are assumed to be read consecutively.
      *
      * \return true if this frame exists, false if at end
      */
      bool readFrame();

      /**
      * Close trajectory file.
      */
      void close();

      private:

          /// Trajectory file
          std::fstream file_;

          /// Number of atoms stored in the file (as declared in file).
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
