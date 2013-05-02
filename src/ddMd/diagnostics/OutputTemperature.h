#ifndef DDMD_OUTPUT_TEMPERATURE_H
#define DDMD_OUTPUT_TEMPERATURE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/diagnostics/Diagnostic.h>
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically write simulation energies to file.
   *
   * \ingroup DdMd_Diagnostic_Module
   */
   class OutputTemperature : public Diagnostic
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      OutputTemperature(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~OutputTemperature()
      {}

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
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
      * Clear nSample counter.
      */
      virtual void clear();

      /**
      * Dump configuration to file
      *
      * \param iStep MD step index
      */
      virtual void sample(long iStep);

   private:

      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long    nSample_;

      /// Has readParam been called?
      long    isInitialized_;
   };

}
#endif
