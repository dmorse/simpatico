#include "SliplinkMcAnalyzerFactory.h"  
#include <mcMd/mcSimulation/McSimulation.h>  
#include <mcMd/mcSimulation/McSystem.h>  

// Include headers for any user defined Analyzers for MC simulations
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
   Analyzer* SliplinkMcAnalyzerFactory::factory(const std::string &className) const
   {
      Analyzer* ptr = 0;

      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "G1MSD") {
         ptr = new G1MSD(system());
      }             
      else if (className == "EndtoEnd") {
         ptr = new EndtoEnd(system());
      }
      else if (className == "EndtoEndXYZ") {
         ptr = new EndtoEndXYZ(system());
      }      
      else if (className == "LinkLengthDist") {
         ptr = new LinkLengthDist(system());
      }  
      else if (className == "LinkLifeTime") {
         ptr = new LinkLifeTime(system());
      }         
      else if (className == "SSChainDist") {
         ptr = new SSChainDist(system());
      }
      else if (className == "VelProf") {
         ptr = new VelProf(system());
      }         
      else if (className == "NLinkAverage") {
         ptr = new NLinkAverage(system());
      }
      else if (className == "InterIntraLink") {
         ptr = new InterIntraLink(system());
      } 
      else if (className == "LinkLTPos") {
         ptr = new LinkLTPos(system());
      }    
      else if (className == "LinkMSD") {
         ptr = new LinkMSD(system());
      }        

      #if 0
      // If not a user-defined class, try the standard factory 
      if (!ptr) {
         ptr = McAnalyzerFactory::factory(className);
      }
      #endif

      return ptr;
   }

}
