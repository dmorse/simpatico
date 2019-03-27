#ifndef MCMD_MC_ENERGY_ANALYZER_H
#define MCMD_MC_ENERGY_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/AverageListAnalyzer.h>
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
   * total simulation energy (kinetic, pair, bond, etc.) as well as
   * for the total, and periodically outputs block averages of each
   * to a file.
   *
   * \sa \ref mcMd_analyzer_McEnergyAnalyzer_page "param file format"
   *
   * \ingroup McMd_MdAnalyzer_Module
   */
   class McEnergyAnalyzer : public AverageListAnalyzer<McSystem>
   {

   public:
   
      /**
      * Constructor.
      *
      * \param system  parent McSystem object. 
      */
      McEnergyAnalyzer(McSystem& system);
   
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);

   protected:

      /**
      * Compute values of energy components, store in values_ array.
      */  
      void compute();
 
   private: 
 
      #ifndef SIMP_NOPAIR
      /// Array index for pair energy accumulator.
      int pairId_;
      #endif

      #ifdef SIMP_BOND
      /// Array index for bond energy accumulator.
      int bondId_;
      #endif

      #ifdef SIMP_ANGLE
      /// Array index for angle energy accumulator.
      int angleId_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Array index for dihedral energy accumulator.
      int dihedralId_;
      #endif

      #ifdef SIMP_COULOMB
      /// Array index for coulomb energy accumulator.
      int coulombId_;
      #endif

      #ifdef SIMP_EXTERNAL
      /// Array index for external energy accumulator.
      int externalId_;
      #endif

      /// Array index for total potential energy accumulator.
      int totalId_;

   };

}
#endif 
