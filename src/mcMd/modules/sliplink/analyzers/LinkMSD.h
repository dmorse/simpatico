#ifndef LINK_MSD_H
#define LINK_MSD_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/simulation/System.h>             // class template parameter
#include <util/accumulators/Average.h>          // member
#include <util/containers/DArray.h>             // member template
#include <util/space/Vector.h>                   // member template parameter
#include <mcMd/links/LinkMaster.h>
#include <mcMd/links/LinkEvents.h>
#include <util/misc/Observer.h>
#include <util/global.h>

#include <cstdio> 
#include <cstring> 

namespace McMd {
   class Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   *
   * Evaluates msd of link ends along the chains
   *
   * \ingroup Analyzer_Module
   */
   class LinkMSD : public SystemAnalyzer<System>,
                   public Observer<LinkAddEvent>, 
                   public Observer<LinkRemoveEvent>  
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      LinkMSD(System &system);
   
      /**
      * Read parameters from file, and allocate data arrays.
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Empty.
      */
      virtual void setup();
   
      /**
      * Evaluate.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);
      
      virtual void update(const LinkAddEvent& event);
      
      virtual void update(const LinkRemoveEvent& event);      
   
      /**
      * Output param file at end of simulation.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Size of the accumulator array
      int capacity_;
      
      /// Total number of links allowed in the system
      int linkCapacity_;
      
      /// Number of samples per block average output.
      int nSamplePerBlock_;
      
      /// Array of Average objects - statistical accumulators
      DArray<Average>  accumulator_;
            
      /// i0_[i] contains the position of each link end at t0
      DArray<int> i0_;
      
      /// t0_[i] contains the initial time for every link end during its life
      DArray<int> t0_;    
         
   };

}
#endif
