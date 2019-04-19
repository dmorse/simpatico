#ifndef MCMD_MD_PRESSURE_ANALYZER_H
#define MCMD_MD_PRESSURE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/AverageListAnalyzer.h>
#include <mcMd/mdSimulation/MdSystem.h>

namespace Util{
   class Average;
}

namespace McMd
{

   using namespace Util;

   /**
   * Compute averages and output block averages of pressure components.
   *
   * This class computes separate averages for each component of the
   * total simulation pressure (kinetic, pair, bond, etc.) as well as
   * for the total, and periodically outputs block averages of each
   * to a file.
   *
   * \sa \ref mcMd_analyzer_MdPressureAnalyzer_page "param file format"
   *
   * \ingroup McMd_MdAnalyzer_Module
   */
   class MdPressureAnalyzer : public AverageListAnalyzer<MdSystem>
   {

   public:
   
      /**
      * Constructor.
      *
      * \param system  parent MdSystem object. 
      */
      MdPressureAnalyzer(MdSystem& system);
   
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load parameters from archive when restarting. 
      *
      * \param ar loading/input archive
      */
      virtual void loadParameters(Serializable::IArchive& ar); 
   
   protected:

      /**
      * Compute and store values of pressure components.
      */  
      void compute();
 
   private: 
 
      /// Array index for virial pressure accumulator.
      int virialId_;

      /// Array index for kinetic pressure accumulator.
      int kineticId_;

      /// Array index for total pressure accumulator.
      int totalId_;

   };

}
#endif 
