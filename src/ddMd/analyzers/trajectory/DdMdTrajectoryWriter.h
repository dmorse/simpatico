#ifndef DDMD_DDMD_TRAJECTORY_WRITER_H
#define DDMD_DDMD_TRAJECTORY_WRITER_H

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
      * Write trajectory file header.
      *
      * \param file output file stream
      */
      void writeHeader(std::ofstream &file);

      /**
      * Write a single frame. 
      *
      * \param file output file stream
      * \param iStep MD time step index
      */
      void writeFrame(std::ofstream &file, long iStep);

   private:

      /// Number of atoms in the file.
      int nAtom_;

      #if 0
      /**
      * Private method to save Group<N> objects.
      */
      template <int N>
      int writeGroups(BinaryFileOArchive& ar,
                      GroupStorage<N>& storage, 
                      GroupCollector<N>& collector); 
      #endif

   };

}
#endif
