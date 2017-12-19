#ifndef LINK_LT_POS_CPP
#define LINK_LT_POS_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinkLTPos.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>
#include <simp/boundary/Boundary.h>
#include <util/misc/FileMaster.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /// Constructor.
   LinkLTPos::LinkLTPos(System& system) 
    : SystemAnalyzer<System>(system)
   {  setClassName("LinkLTPos"); }

   /// Read parameters from file.
   void LinkLTPos::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "linkCapacity", linkCapacity_);      
      readBegin(in,"Distribution");
      read<double>(in, "min", min_);
      read<double>(in, "max", max_);
      read<int>(in,   "nBin", nBin_);
      readEnd(in);
           
   }

   
   /// Register Observer.
   void LinkLTPos::setup() 
   {  
     system().linkMaster().Notifier<LinkAddEvent>::registerObserver(*this);
     system().linkMaster().Notifier<LinkRemoveEvent>::registerObserver(*this); 
     system().linkMaster().Notifier<ReSetAtomEvent>::registerObserver(*this);

     int i, j, chLength = system().simulation().species(0).nAtom();    
     
     birthTimes_.allocate(2*linkCapacity_, chLength);
     
     //Allocate acumulator_ array.
     accumulator_.allocate(chLength);     
  
     for (i=0; i<chLength; ++i) {
        accumulator_[i].setParam(min_, max_, nBin_);  
     } 
          
     for (j=0; j<2*linkCapacity_; ++j) {
       for (i=0; i<chLength; ++i) {
         birthTimes_(j,i)=-1;
       }
     }     
     
   }
 
   /// Don't do anything.
   void LinkLTPos::sample(long iStep) 
   { }  

   void LinkLTPos::update(const LinkAddEvent& event)
   { 
     birthTimes_((event.get()->tag())*2,event.get()->atom0().indexInMolecule()) = system().simulation().iStep();
     birthTimes_((event.get()->tag())*2+1,event.get()->atom1().indexInMolecule()) = system().simulation().iStep();     
   }

   void LinkLTPos::update(const LinkRemoveEvent& event)
   { 
     int btime0, btime1, dtime, ltime;
     int i,chLength = system().simulation().species(0).nAtom();
     int le, le1;
     
     dtime = system().simulation().iStep();     
     le = (event.get()->tag())*2;
     le1 = le+1;
     
     for (i=0; i<chLength; ++i){
       btime0 = birthTimes_(le,i);
       btime1 = birthTimes_(le1,i);            
       if (btime0!=-1) {
         ltime = dtime - btime0;
         birthTimes_(le,i) = -1; 
         accumulator_[i].sample((double)ltime); 
       }
       if (btime1!=-1) {
         ltime = dtime - btime1;
         birthTimes_(le1,i) = -1;
         accumulator_[i].sample((double)ltime); 
       }                  
     }
   }
   
     void LinkLTPos::update(const ReSetAtomEvent& event)
   { 
     
     int tag = event.getLink()->tag();
     int idAtom0 = event.getLink()->atom0().indexInMolecule();
     int idAtom1 = event.getLink()->atom1().indexInMolecule();
     int end = event.getEndId(); 
     int le;
     
     if (end==0) {
       le = tag*2;
       if (birthTimes_(le,idAtom0)==-1) { 
         birthTimes_(le,idAtom0) = system().simulation().iStep();
       }
     }
     if (end==1) {
       le = tag*2+1;
       if (birthTimes_(le,idAtom1)==-1) {
         birthTimes_(le,idAtom1) = system().simulation().iStep();
       }
     }
     
   }

   /// Output results to file after simulation is completed.
   void LinkLTPos::output() 
   {  

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output data to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      int i, chLength = system().simulation().species(0).nAtom();
      for (i=0; i<chLength; ++i) {
        accumulator_[i].output(outputFile_);
        outputFile_ << std::endl;
      }
      outputFile_.close();

   }

}
#endif
