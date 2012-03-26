#ifndef MCMD_VOLUME_AVERAGE_H
#define MCMD_VOLUME_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <util/accumulators/Average.h>         // member template 
#include <mcMd/util/FileMaster.h>
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;

   /**
   * Autocorrelation for vector separation of two atoms on a molecule.  
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class VolumeAverage : public SystemDiagnostic<System>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      VolumeAverage(System &system);
  
      /** 
      * Read parameters from file.
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);
  
      /** 
      * Clear accumulator.
      */
      virtual void initialize();
   
      /** 
      * Evaluate volume of simulation cell, and add to ensemble.
      *
      * \param iStep counter for number of steps
      */
      virtual void sample(long iStep);
   

      /**
      * Output results to file after simulation is completed.
      */
      virtual void output();

      /**
      * Save state to binary file archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchiveType& ar);

      /**
      * Load state from a binary file archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void load(Serializable::IArchiveType& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Array of Average objects - statistical accumulators
      DArray<Average>  accumulators_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Has readParam been called?
      bool    isInitialized_;
   
   };
   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void VolumeAverage::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }
      for (int i = 0; i < Dimension+1; ++i) {
         ar & accumulators_[i];
      }
   }

}
#endif
