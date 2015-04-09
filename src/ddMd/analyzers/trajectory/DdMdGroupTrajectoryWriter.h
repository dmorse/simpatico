#ifndef DDMD_DDMD_GROUP_TRAJECTORY_WRITER_H
#define DDMD_DDMD_GROUP_TRAJECTORY_WRITER_H

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
   * Native binary trajectory format for ddSim, for one group.
   *
   * \ingroup McMd_TrajectoryWriter_Module
   */
   class DdMdGroupTrajectoryWriter : public TrajectoryWriter
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object
      */
      DdMdGroupTrajectoryWriter(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~DdMdGroupTrajectoryWriter();

      /**
      * Read parameter file block.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Read trajectory file header and initialize simulation parameters.
      *
      * \param file output file stream
      */
      void writeHeader(std::ofstream &file);

      /**
      * Read a single frame. Frames are assumed to be read consecutively.
      *
      * \param file output file stream
      * \param iStep MD time step index
      */
      void writeFrame(std::ofstream &file, long iStep);

   private:

      /// Integer id of atom group.
      unsigned int groupId_;

      /// Number of atoms in the file.
      int nAtom_;

   };

}
#endif
