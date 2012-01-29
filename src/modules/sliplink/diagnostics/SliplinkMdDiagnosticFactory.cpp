#ifndef SLIPLINK_MD_DIAGNOSTIC_FACTORY_CPP
#define SLIPLINK_MD_DIAGNOSTIC_FACTORY_CPP

#include "SliplinkMdDiagnosticFactory.h"       // Class header
#include <mcMd/mdSimulation/MdSystem.h>  

// Include headers for any user defined Diagnostics for MD simulations
#include "Crosslinker.h"
#include "G1MSD.h"
#include "EndtoEnd.h"
#include "EndtoEndXYZ.h"
#include "LinkLengthDist.h"
#include "LinkLifeTime.h"
#include "SSChainDist.h"
#include "VelProf.h"
#include "NLinkAverage.h"
#include "InterIntraLink.h"
#include "LinkLTPos.h"
#include "LinkMSD.h"

namespace McMd
{

   /* 
   * Return a pointer to a new instance of className.
   */
   Diagnostic* SliplinkMdDiagnosticFactory::factory(const std::string &className) const
   {
      Diagnostic* spp = 0;

      // Check names of user defined subclasses of Diagnostic
      //if (className == "NewDiagnostic1") {
      //   spp = new NewDiagnostic1(system());
      //} else 
      //if (className == "NewDiagnostic2") {
      //   spp = new NewDiagnostic2(system());
      //} else 
      //  ...
      //}
      
      if (className == "Crosslinker") {
         spp = new Crosslinker(system());
      } 
      else if (className == "G1MSD") {
         spp = new G1MSD(system());
      }  
      else if (className == "EndtoEnd") {
         spp = new EndtoEnd(system());
      }
      else if (className == "EndtoEndXYZ") {
         spp = new EndtoEndXYZ(system());
      }           
      else if (className == "LinkLengthDist") {
         spp = new LinkLengthDist(system());
      } 
      else if (className == "LinkLifeTime") {
         spp = new LinkLifeTime(system());
      }               
      else if (className == "SSChainDist") {
         spp = new SSChainDist(system());
      }  
      else if (className == "VelProf") {
         spp = new VelProf(system());
      }     
      else if (className == "NLinkAverage") {
         spp = new NLinkAverage(system());
      }
      else if (className == "InterIntraLink") {
         spp = new InterIntraLink(system());
      }  
      else if (className == "LinkLTPos") {
         spp = new LinkLTPos(system());
      }    
      else if (className == "LinkMSD") {
         spp = new LinkMSD(system());
      }                         
     
      // If not a user-defined class, try the standard factory 
      if (!spp) {
         spp = MdDiagnosticFactory::factory(className);
      }

      return spp;
   }

}

#endif
