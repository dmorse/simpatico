#ifndef MCMD_MC_WRITE_RESTART_H
#define MCMD_MC_WRITE_RESTART_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/Diagnostic.h>
#include <mcMd/mcSimulation/McSimulation.h>

namespace McMd
{

   using namespace Util;

   class McSimulation;

   /**
   * Periodically write restart file.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class McWriteRestart : public Diagnostic
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent McSimulation object. 
      */
      McWriteRestart(McSimulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~McWriteRestart()
      {} 
   
      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      // Use inherited Diagnostic::loadParameters and Diagnostic::save
   
      /**
      * Write restart file.
      *
      * \param iStep MC step index
      */
      void sample(long iStep);
  
   private:
      
      // Pointer to parent McSimulation
      McSimulation* simulationPtr_;

   };

}
#endif 
