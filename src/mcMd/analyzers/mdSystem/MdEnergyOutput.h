#ifndef MCMD_MD_ENERGY_OUTPUT_H
#define MCMD_MD_ENERGY_OUTPUT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/mdSimulation/MdSystem.h>         // base template parameter
#include <util/global.h> 

namespace McMd
{

   using namespace Util;

   /**
   * Analyzer to output total potential and kinetic energies.
   * 
   * \ingroup McMd_Analyzer_Module
   */
   class MdEnergyOutput : public SystemAnalyzer<MdSystem>
   {

   public:
  
      /* 
      * Constructor.
      */
      MdEnergyOutput(MdSystem& system);

      /**
      * Read output file and nStepPerSample.
      *
      * \param in input parameter stream
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
      * Evaluate energy and print.
      */
      virtual void sample(long iStep);

      /**
      * Output final summary and file format.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;
   
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void MdEnergyOutput::serialize(Archive& ar, const unsigned int version)
   {  Analyzer::serialize(ar, version); }

}
#endif 
