#ifndef LINK_LT_POS_H
#define LINK_LT_POS_H

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
   * LinkLTPos evaluates how long the slip-springs live as a function of position
   * along the chain.
   *
   * \ingroup Analyzer_Module
   */
   class LinkLTPos : public SystemAnalyzer<System>,
                        public Observer<LinkAddEvent>, 
                        public Observer<LinkRemoveEvent>,
                        public Observer<ReSetAtomEvent>
   {
   
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      LinkLTPos(System &system);

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
      
       virtual void update(const ReSetAtomEvent& event);     

      /** 
      * Output results to output file.
      */
      virtual void output();

   private:

      // Output file stream
      std::ofstream outputFile_;
      
      // Distribution statistical accumulator
      DArray<Distribution>  accumulator_;
      
      // Parameters for the statistical accumulator;
      double  min_;
      double  max_;
      int     nBin_;
      
      // Dinamic 2d-array of creation times for each link end le and position
      // i along the chain le=tag*2+endID 0<=le<2*linkCapacity, 0<=i<nMolecule
      DMatrix<int>  birthTimes_;
      
      // Allocated dimension for birthTimes_ array, same as links_ container.
      int linkCapacity_;
  
   };

}
#endif
