#ifndef MCMD_CONFIG_WRITER_H
#define MCMD_CONFIG_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>
#include <mcMd/simulation/System.h>

namespace McMd
{

   using namespace Util;

   /**
   * Periodically write snapshots to a lammps dump (i.e., trajectory) file
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class LammpsDumpWriter : public SystemAnalyzer<System>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent System object. 
      */
      LammpsDumpWriter(System& system);
   
      /**
      * Destructor.
      */
      virtual ~LammpsDumpWriter()
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
      * Write a frame/snapshot to file
      *
      * \param iStep MC step index
      */
      virtual void sample(long iStep);
 
      /**
      * Close file.
      */ 
      virtual void output(); 

   private:
      
      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;
   
      /// Has readParam been called?
      long isInitialized_;
   
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void LammpsDumpWriter::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nSample_;
   }

}
#endif 
