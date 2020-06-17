#ifndef DDMD_DDMD_GROUP_TRAJECTORY_WRITER_H
#define DDMD_DDMD_GROUP_TRAJECTORY_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
      * Load internal state from an input archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an output archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Write trajectory file header.
      */
      void writeHeader();

      /**
      * Write a single frame. 
      *
      * \param iStep MD time step index
      */
      void writeFrame(long iStep);

   private:

      /// Integer id of atom group.
      unsigned int groupId_;

      /// Number of atoms in the file.
      int nAtom_;

   };

}
#endif
