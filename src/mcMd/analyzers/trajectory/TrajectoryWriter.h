#ifndef MCMD_CONFIG_WRITER_H
#define MCMD_CONFIG_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>
#include <mcMd/simulation/System.h>

namespace McMd
{

   using namespace Util;

   /**
   * Periodically write snapshots to a lammps dump (i.e., trajectory) file
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class TrajectoryWriter : public SystemAnalyzer<System>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent System object. 
      */
      TrajectoryWriter(System& system);
   
      /**
      * Destructor.
      */
      virtual ~TrajectoryWriter()
      {} 
   
      /**
      * Read interval.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Clear nSample counter.
      */
      virtual void setup();
  
      /**
      * Write a frame/snapshot to trajectory file.
      *
      * \param iStep step index
      */
      virtual void sample(long iStep);
  
      /**
      * Close trajectory file after run.
      */
      virtual void output();
  
   protected:
      
      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;
   
      /// Has readParam been called?
      long isInitialized_;

   protected:

      /**
      * Is the file format binary (true) or text (false)?
      */
      bool isBinary() const;

      /**
      * Write data that should appear once, at beginning of the file. 
      *
      * Called by sample on first invocation. Default implementation is empty.
      *
      * \param out output file stream
      */
      virtual void writeHeader(std::ofstream& out)
      {};

      /**
      * Write data that should appear in every frame.
      * 
      * \param out output file stream
      * \param iStep MD time step index
      */
      virtual void writeFrame(std::ofstream& out, long iStep) = 0;

      /// Is the trajectory file a binary file? 
      bool isBinary_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void TrajectoryWriter::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nSample_;
   }

}
#endif 
