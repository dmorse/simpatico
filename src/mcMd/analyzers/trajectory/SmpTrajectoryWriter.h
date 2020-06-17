#ifndef MCMD_SMP_TRAJECTORY_WRITER_H
#define MCMD_SMP_TRAJECTORY_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"

namespace McMd
{

   using namespace Util;

   /**
   * Write a DdMd trajectory file. 
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class SmpTrajectoryWriter : public TrajectoryWriter
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent MdSystem object.
      */
      SmpTrajectoryWriter(MdSystem& system);

      /**
      * Destructor.
      */
      virtual ~SmpTrajectoryWriter();

      /**
      * Serialize to/from an archive.
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   protected:

      /**
      * Open file and write header.
      */
      virtual void writeHeader();

      /**
      * Write data that should appear in every frame.
      *
      * \param iStep MD time step index
      */
      virtual void writeFrame(long iStep);

   private:

      int nAtom_;

   };

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void SmpTrajectoryWriter::serialize(Archive& ar, 
                                        const unsigned int version)
   {  TrajectoryWriter::serialize(ar, version); }

}
#endif
