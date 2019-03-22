#ifndef MCMD_DISTRIBUTION_ANALYZER_H
#define MCMD_DISTRIBUTION_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>    // base class template
#include <util/accumulators/Distribution.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * DistributionAnalyzer evaluates the distribution of a real variable. 
   *
   * \sa \ref mcMd_analyzer_DistributionAnalyzer_page "parameter file format"
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   template <class SystemType>
   class DistributionAnalyzer : public SystemAnalyzer<SystemType>
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent SystemType object
      */
      DistributionAnalyzer(SystemType &system);

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
      * Output results to output file.
      */
      virtual void output();

      using ParamComposite::read;
      using ParamComposite::loadParameter;
      using Analyzer::writeParam;
      using Analyzer::outputFileName;
      using Analyzer::fileMaster;

   protected:

      /// Minimum of range
      double min_;

      /// Maximum of range
      double max_;

      /// Number of bins in range
      double nBin_;

      /** 
      * Add a value to the accumulator.
      */
      void increment(double value);

      /**
      * Get Distribution accumulator by const reference.
      */
      const Distribution& accumulator() const;

      using ParamComposite::setClassName;
      using Analyzer::readInterval;
      using Analyzer::readOutputFileName;
      using Analyzer::loadInterval;
      using Analyzer::loadOutputFileName;

   private:

      // Output file stream
      std::ofstream outputFile_;

      // Distribution statistical accumulator
      Distribution  accumulator_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class SystemType>
   template <class Archive>
   void DistributionAnalyzer<SystemType>::serialize(Archive& ar, 
                                                    const unsigned int version)
   {   
      Analyzer::serialize(ar, version);
      ar & min_;
      ar & max_;
      ar & nBin_;
      ar & accumulator_;
   }

   /*
   * Get Distribution accumulator by const reference.
   */
   template <class SystemType>
   inline
   const Distribution& DistributionAnalyzer<SystemType>::accumulator() const
   {  return accumulator_; }

}
#endif
