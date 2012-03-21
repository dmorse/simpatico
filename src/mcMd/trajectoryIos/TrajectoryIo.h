#ifndef MCMD_TRAJECTORY_IO_H
#define MCMD_TRAJECTORY_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/boundary/Boundary.h>
#include <util/global.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   class Simulation;
   class System;
   
   /**
   * Trajectory file reader and writer. A subclass of TrajectoryIo needs to implement a
   * the pair readHeader() and readFrame() or the pair writeHedaer() and writeFrame() or both.
   * If either pair is not implemented, calling the standard implementation will result in an
   * error message.
   *
   * \ingroup TrajectoryIo_Module
   */
   class TrajectoryIo
   {
   
   public:

      /// Constructor. 
      TrajectoryIo(System& system);
 
      /// Destructor.   
      virtual ~TrajectoryIo();
 
      /**
      * Read trajectory file header. This is called once to initialize the file reader.
      *
      * \param file input file stream.
      */
      virtual void readHeader(std::fstream &file);
 
      /**
      * Read a single frame. Frames are assumed to be read consecutively. Before a call to
      * readFrame(), the readHeader() method must be called.
      *
      * \param file input file stream
      */
      virtual void readFrame(std::fstream& file);

      /**
      * Write trajectory file header.
      *
      * \param file output file stream.
      */
      virtual void writeHeader(std::fstream &file);
 
      /**
      * Write a single frame. Frames are assumed to be written consecutively.
      *
      * \param file output file stream
      */
      virtual void writeFrame(std::fstream& file);

      /**
      * Return number of frames in trajectory file.
      */
      int nFrames() const;
 
   protected:

      /// The number of frames in the trajectory file
      int nFrames_;

      /// Get a reference to the parent System. 
      System &system() const;

      /// Get a reference to the parent Simulation. 
      Simulation &simulation() const;

      /// Get the Boundary.
      Boundary &boundary() const;

   private:
   
      /// Boundary object.
      Boundary   *boundaryPtr_;
   
      /// Pointer to parent System;
      System     *systemPtr_;
   
      /// Pointer to parent Simulation.
      Simulation *simulationPtr_;
   
   }; // end class TrajectoryIo


   // Inline functions 

   /* 
   * Get the number of frames
   */
   inline int TrajectoryIo::nFrames() const
   {
      return nFrames_;
   }

   /* 
   * Get the parent System. 
   */
   inline System& TrajectoryIo::system() const
   {
      assert(systemPtr_); 
      return *systemPtr_; 
   }
 
   /* 
   * Get the parent Simulation.
   */
   inline Simulation& TrajectoryIo::simulation() const
   { 
      assert(simulationPtr_);
      return *simulationPtr_; 
   }

   /* 
   * Get the Boundary.
   */
   inline Boundary& TrajectoryIo::boundary() const
   {
      assert(boundaryPtr_); 
      return *boundaryPtr_; 
   }


} 
#endif
