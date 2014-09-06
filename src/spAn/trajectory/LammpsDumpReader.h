#ifndef SPAN_LAMMPS_DUMP_READER_H
#define SPAN_LAMMPS_DUMP_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/trajectory/TrajectoryReader.h>  // base class

namespace SpAn
{

   class Configuration;
   using namespace Util;

   /**
   * Reader for lammps dump trajectory file format.
   *
   * \ingroup SpAn_Trajectory_Module
   */
   class LammpsDumpReader  : public TrajectoryReader
   {

   public:

      /**
      * Default constructor.
      */
      LammpsDumpReader();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object
      */
      LammpsDumpReader(Configuration& configuration);

      /**
      * Destructor.
      */
      virtual ~LammpsDumpReader();

      /**
      * Read a frame.
      *
      * \param file input file 
      * \return true if a frame was found, false if end of file
      */
      virtual bool readFrame(std::ifstream& file);

   };

}
#endif
