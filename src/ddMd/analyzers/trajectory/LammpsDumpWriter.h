#ifndef DDMD_LAMMPS_DUMP_WRITER_H
#define DDMD_LAMMPS_DUMP_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/trajectory/TrajectoryWriter.h>   // base class

namespace DdMd
{

   using namespace Util;

   /**
   * Write a trajectory in the Lammps dump format.
   *
   * \ingroup DdMd_Analyzer_Trajectory_Module
   */
   class LammpsDumpWriter : public TrajectoryWriter
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object
      */
      LammpsDumpWriter(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~LammpsDumpWriter();

      /**
      * Read a single frame. Frames are assumed to be read consecutively.
      *
      * \param file output file stream
      * \param iStep MD time step index
      */
      void writeFrame(std::ofstream &file, long iStep);

   private:

      /// Number of atoms in the file.
      int nAtom_;

   };

}
#endif
