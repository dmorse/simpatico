#ifndef MCMD_TYPE_DISTRIBUTION_H
#define MCMD_TYPE_DISTRIBUTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>      // base class template
#include <mcMd/mcSimulation/McSystem.h>         // base template parameter
#include <util/accumulators/IntDistribution.h>  // member

namespace Simp {
   class Species;
}

namespace McMd
{

   class SpeciesMutator;

   using namespace Util;
   using namespace Simp;

   /**
   * Calculate distribution of type indices for mutable species.
   *
   * \ingroup McMd_Analyzer_Mc_Module
   */
   class TypeDistribution : public SystemAnalyzer<McSystem>
   {

   public:
  
      /**
      * Constructor.
      */
      TypeDistribution(McSystem& system);

      /**
      * Read output file and nStepPerSample.
      *
      * \param in input parameters stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /* 
      * Evaluate energy and print.
      */
      void sample(long iStep);

      /*
      * Output final summary and file format.
      */
      virtual void output();

   private:

     /// Output file stream.
     std::ofstream outputFile_;

     /// Histogram of values.
     DArray<int> distribution_;

     /// Number of samples.
     int nSample_;

     /// Species index.
     int speciesId_;

     /// Number of possible internal states.
     int nState_;

     /// Pointer to Species.
     Species* speciesPtr_;

     /// Pointer to associated SpeciesMutator.
     SpeciesMutator* mutatorPtr_;
   
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void TypeDistribution::serialize(Archive& ar, const unsigned int version)
   {  
       Analyzer::serialize(ar, version); 
       ar & speciesId_;
       ar & nState_;
       ar & nSample_;
       ar & distribution_;
   }

}
#endif 
