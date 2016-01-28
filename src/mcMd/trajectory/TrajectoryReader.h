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
      * Note: By convention, this function does not add any prefixes
      * to the trajectory file path. The filename argument should thus 
      * be given as a relative path defined relative to the directory 
      * from which the program is executed.
      *
      * \param filename input file name, relative to working directory.
      */
      virtual void open(std::string filename) = 0;

      /**
      * Read a single frame. Frames are assumed to be read consecutively. 
      *
      * \return true if a frame is avaiable, false if at end of file
      */
      virtual bool readFrame() = 0;

      /**
      * Close trajectory file.
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
