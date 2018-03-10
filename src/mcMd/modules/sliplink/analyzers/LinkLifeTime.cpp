#ifndef LINK_LIFE_TIME_CPP
#define LINK_LIFE_TIME_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinkLifeTime.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>
#include <simp/species/Species.h>
#include <simp/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/misc/FileMaster.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /// Constructor.
   LinkLifeTime::LinkLifeTime(System& system) 
    : SystemAnalyzer<System>(system)
   {  setClassName("LinkLifeTime"); }

   /// Read parameters from file.
   void LinkLifeTime::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "linkCapacity", linkCapacity_);      
      readParamComposite(in, accumulator_);
   }

   
   /// Register Observer.
   void LinkLifeTime::setup() 
   {  
     system().linkMaster().Notifier<LinkAddEvent>::registerObserver(*this);
     system().linkMaster().Notifier<LinkRemoveEvent>::registerObserver(*this); 
     birthTimes_.allocate(linkCapacity_);
     accumulator_.clear();
     for (int i=0; i<linkCapacity_; ++i) {
       birthTimes_[i]=-1;
     }
   }
 
   /// Don't do anything.
   void LinkLifeTime::sample(long iStep) 
   { }  

   void LinkLifeTime::update(const LinkAddEvent& event)
   { 
     birthTimes_[event.get()->tag()] = system().simulation().iStep();
   }

   void LinkLifeTime::update(const LinkRemoveEvent& event)
   { 
     int btime, dtime, ltime;
     
     btime = birthTimes_[event.get()->tag()];
     dtime = system().simulation().iStep();
                      
     if (btime!=-1) {
       ltime = dtime - btime;
       birthTimes_[event.get()->tag()] = -1;
       accumulator_.sample((double)ltime); 
     }
     
   }

   /// Output results to file after simulation is completed.
   void LinkLifeTime::output() 
   {  

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output data to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
#endif
