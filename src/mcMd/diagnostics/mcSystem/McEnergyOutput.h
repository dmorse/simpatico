#ifndef MCMD_MC_ENERGY_OUTPUT_H
#define MCMD_MC_ENERGY_OUTPUT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h> // base class template
#include <mcMd/mcSimulation/McSystem.h>        // base template parameter

namespace McMd
{

   using namespace Util;

   /**
   * Diagnostic to output total potential energy.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class McEnergyOutput : public SystemDiagnostic<McSystem>
   {

   public:

      /**   
      * Constructor.
      *
      * \param system parent McSystem.
      */
      McEnergyOutput(McSystem& system);

      /**
      * Read output file and nStepPerSample.
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
      * Output final summary and file format
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
   void McEnergyOutput::serialize(Archive& ar, const unsigned int version)
   {  Diagnostic::serialize(ar, version); }

}
#endif 
