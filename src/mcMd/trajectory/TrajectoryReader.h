#ifndef MCMD_TRAJECTORY_READER_H
#define MCMD_TRAJECTORY_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/boundary/Boundary.h>  // typedef
#include <util/global.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   class Simulation;
   class System;

   /**
   * Trajectory file reader (base class).
   *
   * \ingroup McMd_Trajectory_Module
   */
   class TrajectoryReader
   {

   public:

      /**
      * Constructor.
      */
      TrajectoryReader(System& system);

      /**
      * Destructor.
      */
      virtual ~TrajectoryReader();

      /**
      * Open trajectory file and read header, if any.
      *
      * By convention, this function treats the trajectory filename
      * as the name of an input file, and opens the file using the 
      * FileMaster:openInutFile function. This function prepends the 
      * input prefix (if any) to the file path. If compiled with MPI 
      * enabled, so that each processor simulates a different system, 
      * it also prepends a processor id prefix before the input prefix.
      *
      * \param filename trajectory input file name.
      */
      virtual void open(std::string filename) = 0;

      /**
      * Read a single frame. Frames are assumed to be read consecutively. 
      *
      * This function reads a frame from the trajectory file that was
      * opened by the open() function.
      *
      * \return true if a frame is avaiable, false if at end of file
      */
      virtual bool readFrame() = 0;

      /**
      * Close the trajectory file.
      */
      virtual void close() = 0;

   protected:

      /// Total number of atoms (all species)
      int nAtomTotal_;

      /// Get a reference to the parent System.
      System &system() const;

      /// Get a reference to the parent Simulation.
      Simulation &simulation() const;

      /// Get the Boundary.
      Boundary &boundary() const;

      /**
      * Add all molecules to system. 
      */
      virtual void addMolecules();

   private:

      /// Boundary object.
      Boundary *boundaryPtr_;

      /// Pointer to parent System;
      System *systemPtr_;

      /// Pointer to parent Simulation.
      Simulation *simulationPtr_;

   }; // end class TrajectoryReader

   // Inline functions

   /*
   * Get the parent System.
   */
   inline System& TrajectoryReader::system() const
   {
      assert(systemPtr_);
      return *systemPtr_;
   }

   /*
   * Get the parent Simulation.
   */
   inline Simulation& TrajectoryReader::simulation() const
   {
      assert(simulationPtr_);
      return *simulationPtr_;
   }

   /*
   * Get the Boundary.
   */
   inline Boundary& TrajectoryReader::boundary() const
   {
      assert(boundaryPtr_);
      return *boundaryPtr_;
   }

}
#endif
