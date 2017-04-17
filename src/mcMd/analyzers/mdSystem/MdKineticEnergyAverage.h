#ifndef MCMD_MD_KINETIC_ENERGY_AVERAGE_H
#define MCMD_MD_KINETIC_ENERGY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>   // base class template
#include <mcMd/mdSimulation/MdSystem.h>          // base template parameter
#include <util/accumulators/Average.h>           // member

#include <fstream> 

namespace McMd
{

   using namespace Util;

   /**
   * MdKineticEnergyAverage averages of total kinetic energy.
   *
   * See \ref mcMd_analyzer_MdKineticEnergyAverage_page "here" for 
   * the parameter file format and any other user documentation.
   *
   * \ingroup McMd_Analyzer_Md_Module
   */
   class MdKineticEnergyAverage : public SystemAnalyzer<MdSystem>
   {
   
   public:

      /**   
      * Constructor.
      */
      MdKineticEnergyAverage(MdSystem& system);

      /**
      * Read parameters and initialize.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);
   
      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
   
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
   
      /* 
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      // Has readParam been called?
      bool isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void MdKineticEnergyAverage::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & accumulator_;
   }

}
#endif 
