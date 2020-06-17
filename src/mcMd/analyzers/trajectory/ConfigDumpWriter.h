#ifndef MCMD_LAMMPS_DUMP_WRITER_H
#define MCMD_LAMMPS_DUMP_WRITER_H

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
   * Periodically write snapshots to a lammps dump (i.e., trajectory) file
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class ConfigDumpWriter : public TrajectoryWriter
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object.
      */
      ConfigDumpWriter(System& system);

      /**
      * Destructor.
      */
      virtual ~ConfigDumpWriter();

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

   };

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void ConfigDumpWriter::serialize(Archive& ar, const unsigned int version)
   {  TrajectoryWriter::serialize(ar, version); }

}
#endif
