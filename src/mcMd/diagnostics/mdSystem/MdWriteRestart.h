#ifndef MCMD_MD_WRITE_RESTART_H
#define MCMD_MD_WRITE_RESTART_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/Diagnostic.h>
#include <mcMd/mdSimulation/MdSimulation.h>

namespace McMd
{

   using namespace Util;

   class MdSimulation;

   /**
   * Periodically write restart file.
   *
   * \ingroup Diagnostic_Module
   */
   class MdWriteRestart : public Diagnostic
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent MdSimulation object
      */
      MdWriteRestart(MdSimulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~MdWriteRestart()
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

      // Pointer to parent MdSimulation
      MdSimulation* simulationPtr_;

   };

}
#endif 
