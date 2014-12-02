#ifndef MCMD_BOUNDARY_AVERAGE_H
#define MCMD_BOUNDARY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <util/accumulators/Average.h>          // member template 
#include <util/misc/FileMaster.h>

namespace McMd
{

   using namespace Util;

   /**
   * Average of boundary lengths and volume of simulation cell.  
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class BoundaryAverage : public SystemAnalyzer<System>
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
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);
  
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
      Analyzer::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & accumulators_;
   }

}
#endif
