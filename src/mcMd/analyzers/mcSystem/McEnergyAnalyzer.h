#ifndef MCMD_MC_ENERGY_ANALYZER_H
#define MCMD_MC_ENERGY_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>
#include <mcMd/mcSimulation/McSystem.h>

namespace Util{
   class Average;
}

namespace McMd
{

   using namespace Util;

   /**
   * Compute averages and output block averages of energy components.
   *
   * This class computes separate averages for each component of the
   * total simulation energy (pair, bond, angle, etc.) as well as
   * for the total, and periodically outputs block averages of each
   * to a file.
   *
   * \sa \ref mcMd_analyzer_McEnergyAnalyzer_page "param file format"
   *
   * \ingroup McMd_McAnalyzer_Module
   */
   class McEnergyAnalyzer : public SystemAnalyzer<McSystem>
   {

   public:
   
      /**
      * Constructor.
      *
      * \param system  parent McSystem object. 
      */
      McEnergyAnalyzer(McSystem& system);
   
      /**
      * Destructor.
      */
      virtual ~McEnergyAnalyzer()
      {} 
   
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an input archive.
      *
      * \param ar  input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an output archive.
      *
      * \param ar  output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /**
      * Clear nSample counter and all accumulators.
      */
      virtual void clear();
  
      /**
      * Setup before main loop.
      *
      * Open all files required for output during simulation.
      */
      virtual void setup();

      /**
      * Write energy to file
      *
      * \param iStep MD time step index
      */
      virtual void sample(long iStep);

      /**
      * Write file averages and error analysis to file
      */
      virtual void output();

   private:
 
      // Output file stream
      std::ofstream outputFile_;

      // Pointers to average accumulators for energy and components
      Average* totalAveragePtr_;
      #ifndef SIMP_NOPAIR
      Average* pairAveragePtr_;
      #endif
      #ifdef SIMP_BOND
      Average* bondAveragePtr_;
      #endif
      #ifdef SIMP_ANGLE
      Average* angleAveragePtr_;
      #endif
      #ifdef SIMP_DIHEDRAL
      Average* dihedralAveragePtr_;
      #endif
      #ifdef SIMP_EXTERNAL
      Average* externalAveragePtr_;
      #endif

      // Number of sample per block average
      int nSamplePerBlock_;

      /// Has readParam been called?
      long isInitialized_;
   
   };

}
#endif 
