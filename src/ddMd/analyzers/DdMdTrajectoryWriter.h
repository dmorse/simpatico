#ifndef DDMD_DDMD_TRAJECTORY_WRITER_H
#define DDMD_DDMD_TRAJECTORY_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/TrajectoryWriter.h>   // base class

namespace DdMd
{

   using namespace Util;

   /**
   * Native binary trajectory format for ddSim.
   *
   * \ingroup McMd_TrajectoryWriter_Module
   */
   class DdMdTrajectoryWriter : public TrajectoryWriter
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object
      */
      DdMdTrajectoryWriter(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~DdMdTrajectoryWriter();

      /**
      * Read trajectory file header and initialize simulation parameters.
      *
      * \param file output file stream
      * \param iStep MD time step index
      */
      void writeHeader(std::ofstream &file, long iStep);

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

      /**
      * Private method to save Group<N> objects.
      */
      template <int N>
      int writeGroups(BinaryFileOArchive& ar,
                      GroupStorage<N>& storage, 
                      GroupCollector<N>& collector); 

   };

}
#endif
