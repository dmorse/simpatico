#ifndef MCMD_CLUSTERS_STATISTICS_H
#define MCMD_CLUSTERS_STATISTICS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>              // base class template
#include <mcMd/simulation/System.h>                     // class template parameter
#include <mcMd/neighbor/CellList.h>                     // member
#include <util/accumulators/IntDistribution.h>          // member
#include <util/containers/DArray.h>                     // member template
#include <util/containers/GArray.h>                     // member template
#include <util/containers/FSArray.h>                    // member template
#include <util/space/Vector.h>                          // member template parameter

#include <cstdio> 
#include <cstring> 

namespace McMd
{

   using namespace Util;

   class Species;

   /**
   * This class is intended to identify Clusters in polymeric systems.
   */
   class ClustersStatistics : public SystemAnalyzer<System>
   {
   
   public:
      /**
      * Cluster struct definition.
      */
      struct Cluster
      {
         Molecule* self_;
         int clusterId_;
      };

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      ClustersStatistics(System &system);
   
      /**
      * Read parameters from file, and allocate data array.
      *
      * Input format:
      *
      *   - int    interval        : sampling interval
      *   - string outputFileName  : base name for output file(s)
      *   - int    nSamplePerBlock : interval for output of block averages
      *   - int    speciesId       : integer id for Species of interest
      *
      * No block averages are output if nSamplePerBlock = 0. Otherwise,
      * block averages are output to a file named (outputFileName).dat. 
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /** 
      * Roots out all the Clusters.
      */
      virtual void findClusters(Molecule* molPtr, int clusterId);

      /** 
      * Returns molecule cluster Id
      */
      int clusterId(Molecule* molPtr);

      /**
      * Evaluate squared radii of gyration for all molecules, add to ensemble.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

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
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Clusters species type
      int  speciesId_;

      /// Clusters species type
      Species*  speciesPtr_;

      /// TypeId of touching cores
      int  coreId_;

      /// TypeId of touching cores
      double  cutoff_;

      /// CellList of specified atoms of the species of interest
      CellList cellList_;

      /// Array of relevant species' molecules with Cluster tags
      DArray<Cluster>  clusters_;

      /// Array of length of different Clusters
      GArray<int> clusterLengths_;

      /// Histogram Min bin.
      int  histMin_;

      /// Histogram Min bin.
      int  histMax_;

      /// Distribution of the Clusters.
      IntDistribution  hist_;
   
      /// Number of configurations dumped thus far(first dump is zero).
      long  nSample_;

      /// Has readParam been called?
      bool  isInitialized_;

   };

   /**
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void ClustersStatistics::serialize(Archive& ar, const unsigned int version)
   {  
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & coreId_;
      ar & histMin_;
      ar & histMax_;
   }

}
#endif
