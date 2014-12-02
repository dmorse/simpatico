#ifndef TOOLS_LAMMPS_DUMP_WRITER_H
#define TOOLS_LAMMPS_DUMP_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <tools/analyzers/TrajectoryWriter.h>   // base class

#include <iostream>
#include <fstream>

namespace Tools
{

   class Processor;
   class Configuration;
   using namespace Util;

   /**
   * Write a trajectory in the Lammps dump format.
   *
   * \ingroup Tools_Analyzer_Module
   */
   class LammpsDumpWriter : public TrajectoryWriter
   {

   public:

      /**
      * Constructor.
      *
      * \param processor parent Processor object
      */
      LammpsDumpWriter(Processor& processor);

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object
      * \param fileMaster asssociated Util::FileMaster object 
      */
      LammpsDumpWriter(Configuration& configuration, FileMaster& fileMaster);

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
      void writeFrame(std::ofstream &file, long int iStep);

   private:

      /// Number of atoms in the file.
      int nAtom_;

   };

}
#endif
