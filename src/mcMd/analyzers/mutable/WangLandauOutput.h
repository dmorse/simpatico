#ifndef MCMD_WANG_LANDAU_OUTPUT_H
#define MCMD_WANG_LANDAU_OUTPUT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/mcSimulation/McSystem.h>         // base template parameter
#include <util/accumulators/IntDistribution.h>  // member
#include <util/containers/DArray.h> // member template

namespace McMd
{

   using namespace Util;
   class WangLandauMove;
   /**
   * Save the WangLandau function to file.
   *
   * \ingroup McMd_Analyzer_Mc_Module
   */
   class WangLandauOutput : public SystemAnalyzer<McSystem>
   {

   public:
  
      /** 
      * Constructor.
      */
      WangLandauOutput(McSystem& system);

      /**
      * Read output file and nStepPerSample.
      */
      virtual void readParameters(std::istream& in);
 
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive.
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

      /**
      * Evaluate energy and print.
      */
      void sample(long iStep);

      /**
      * Output final summary and file format.
      */
      virtual void output();

   private:

     /// Output file stream.
     std::ofstream outputFile_;

     /// Histogram of values.
     IntDistribution distribution_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void WangLandauOutput::serialize(Archive& ar, const unsigned int version)
   {  
      Analyzer::serialize(ar, version); 
   }

}
#endif 
