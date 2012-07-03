#ifndef MCMD_BOUNDARY_AVERAGE_H
#define MCMD_BOUNDARY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * Average of boundary lengths and volume of simulation cell.  
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class BoundaryAverage : public SystemDiagnostic<System>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      BoundaryAverage(System &system);
  
      /** 
      * Read parameters from file.
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);
  
      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /** 
      * Evaluate volume and lengths of simulation cell, 
      * and add to ensemble.
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
   
      /// Output file stream.
      std::ofstream outputFile_;

      /// (Dimension + 1) sized array of Average objects. 
      DArray<Average>  accumulators_;

      /// Number of samples per block average object.
      int nSamplePerBlock_;

      /// Has readParam been called?
      bool    isInitialized_;
   
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void BoundaryAverage::serialize(Archive& ar, const unsigned int version)
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
