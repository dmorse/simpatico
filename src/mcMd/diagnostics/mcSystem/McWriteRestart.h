#ifndef MC_WRITE_RESTART_H
#define MC_WRITE_RESTART_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
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
   * \ingroup Diagnostic_Module
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
      * Read interval and file name.
      *
      * \param in input parameter file
      */
      virtual void readParam(std::istream& in);
   
      /**
      * Write restart file.
      *
      * \param iStep MC step index
      */
      void sample(long iStep);
  
   private:
      
      // Name of restart file.
      std::string filename_;

      // Pointer to parent McSimulation
      McSimulation* simulationPtr_;

   };

}
#endif 
