#ifndef DDMD_KINETIC_ENERGY_ANALYZER_H
#define DDMD_KINETIC_ENERGY_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/simulation/Simulation.h>
#include <util/accumulators/Average.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Output and evaluate average of kinetic energy.
   *
   * \sa \ref ddMd_analyzer_KineticEnergyAnalyzer_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class KineticEnergyAnalyzer : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      KineticEnergyAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~KineticEnergyAnalyzer(); 
   
      /**
      * Read dumpPrefix and interval.
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
      * Compute and a sampled value and add it to a sequence.
      *
      * \param iStep MD step index
      */
      virtual void sample(long iStep);

      /**
      * Dump configuration to file
      */
      virtual void output();

   protected:

      /**
      * Function to compute value.
      *
      * Call on all processors.
      */
      virtual void compute();

      /**
      * Current value, set by compute function.
      *
      * Call only on master.
      */
      virtual double& value();

   private:

      /// Output file stream.
      std::ofstream  outputFile_;

      /// Pointer to Average object (only instantiated on master processor)
      Average *accumulatorPtr_;

      /// Current value (set by compute function)
      double value_;
      
      /// Number of samples per block average output.
      int nSamplePerBlock_;
   
      /// Has readParam been called?
      bool isInitialized_;
   
   };

}
#endif 
