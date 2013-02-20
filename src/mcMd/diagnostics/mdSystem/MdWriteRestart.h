#ifndef MCMD_MD_WRITE_RESTART_H
#define MCMD_MD_WRITE_RESTART_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/Diagnostic.h>
#include <mcMd/mdSimulation/MdSimulation.h>

namespace McMd
{

   using namespace Util;

   class MdSimulation;

   /**
   * Repeatedly write restart file.
   *
   * \ingroup McMd_Diagnostic_Module
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
      * Read interval and outputFileName.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      // Using default Diagnostic::loadParameters and Diagnostic::save.
   
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
