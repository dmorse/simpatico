#ifndef LINK_LIFE_TIME_H
#define LINK_LIFE_TIME_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>    // base class template
#include <mcMd/simulation/System.h>               // base class template parameter
#include <util/accumulators/Distribution.h>
#include <util/containers/DArray.h>
#include <mcMd/links/LinkMaster.h>
#include <mcMd/links/LinkEvents.h>
#include <util/misc/Observer.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * LinkLifeTime evaluates how long the slip-springs live.
   *
   * \ingroup Analyzer_Module
   */
   class LinkLifeTime : public SystemAnalyzer<System>,
                        public Observer<LinkAddEvent>, 
                        public Observer<LinkRemoveEvent>
   
   {
   
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      LinkLifeTime(System &system);

      /**
      * Read parameters from file.  
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /** 
      * Register Observer.
      */
      virtual void setup();
   
      /** 
      * Don't do anything.
      *
      * \param iStep step counter
      */
      void sample(long iStep);
      
      virtual void update(const LinkAddEvent& event);
      
      virtual void update(const LinkRemoveEvent& event);

      /** 
      * Output results to output file.
      */
      virtual void output();

   private:

      // Output file stream
      std::ofstream outputFile_;
      
      // Distribution statistical accumulator
      Distribution  accumulator_;
      
      // Dinamic array of creation times for each link
      DArray<int>  birthTimes_;
      
      // Allocated dimension for birthTimes_ array, same as links_ container.
      int linkCapacity_;
  
   };

}
#endif
