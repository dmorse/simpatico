#ifndef LINK_MSD_CPP
#define LINK_MSD_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinkMSD.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>
#include <simp/boundary/Boundary.h>
#include <util/misc/FileMaster.h>        


namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /// Constructor.
   LinkMSD::LinkMSD(System& system) 
    : SystemAnalyzer<System>(system)
   {  setClassName("LinkMSD"); }

   /// Read parameters from file, and allocate data arrays.
   void LinkMSD::readParameters(std::istream& in) 
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"capacity", capacity_);
      read<int>(in,"linkCapacity", linkCapacity_);
      
      //Allocate acumulator_ array.
      accumulator_.allocate(capacity_);
      
      // Allocate arrays 
      int nle = 2*linkCapacity_;
      i0_.allocate(nle);
      t0_.allocate(nle);   
      
      nSamplePerBlock_=0;
      
      int i;
      for (i=0; i < capacity_; ++i){
        accumulator_[i].setNSamplePerBlock(nSamplePerBlock_);
      }
      // If nSamplePerBlock != 0, open an output file for block averages.
      if (nSamplePerBlock_ != 0) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }      
      
   }

   
   /// Empty.
   void LinkMSD::setup() 
   {    
     system().linkMaster().Notifier<LinkAddEvent>::registerObserver(*this);
     system().linkMaster().Notifier<LinkRemoveEvent>::registerObserver(*this); 
     int nle = 2*linkCapacity_;
     int i;
 
     for (i=0; i < capacity_; ++i){
        accumulator_[i].clear();
     }
     for (i=0; i < nle; ++i){
        i0_[i]=-1;
	t0_[i]=-1;
     }
   }
 
   /// Evaluate and output x-velocity profile as a function of z.
   void LinkMSD::sample(long iStep) 
   { 
     if (isAtInterval(iStep))  {
      
       int    nLink = system().linkMaster().nLink();
       int    idLink, dt, dr, dr2, le, le1, i, t;     
       Link*  linkPtr;
               
       for (idLink=0; idLink < nLink; ++idLink){
	 linkPtr = &(system().linkMaster().link(idLink));
	 le = (linkPtr->tag())*2;
	 le1 = le+1;
	 if (t0_[le] != -1) {
	   t = iStep;
	   i = linkPtr->atom0().indexInMolecule();
	   dt = (t - t0_[le])/interval();
	   dr = i - i0_[le];
	   dr2 = dr*dr;
	   //std::cout << "Tag: " << linkPtr->tag() << " end: " << le << " dt: "  << dt << " dr :" << dr << std::endl;
           //std::cout << "i0: " << i0_[le] << " i: " << i << " t0: "  << t0_[le] << " t :" << t << std::endl;	   
	   if (dt < capacity_) {
	     accumulator_[dt].sample(double(dr2), outputFile_);
	   } 
	   /*else {
	      UTIL_THROW("link MSD accumulator capacity exceded");
	   }*/	   
	 }
	 else {
	   t0_[le] = iStep;
	   i0_[le] = linkPtr->atom0().indexInMolecule();
	 }
	 if (t0_[le1] != -1) {
	   t = iStep;
	   i = linkPtr->atom1().indexInMolecule();
	   dt = (t - t0_[le1])/interval();
	   dr = i - i0_[le1];
	   dr2 = dr*dr;
	   //std::cout << "Tag: " << linkPtr->tag() << " end: " << le1 << " dt: "  << dt << " dr :" << dr << std::endl;
           //std::cout << "i0: " << i0_[le1] << " i: " << i << " t0: "  << t0_[le1] << " t :" << t << std::endl;	   
	   if (dt < capacity_) {
	     accumulator_[dt].sample(double(dr2), outputFile_);
	   } 
	   /*else {
	      UTIL_THROW("link MSD accumulator capacity exceded");
	   }*/
	 }
	 else {
	   t0_[le1] = iStep;
	   i0_[le1] = linkPtr->atom1().indexInMolecule();
	 }
       }

     } // if isAtInterval
   }
   
   void LinkMSD::update(const LinkAddEvent& event)
   {
     int le = (event.get()->tag())*2;
     int le1 = le+1;
     //std::cout << "AddEvent Tag: " << event.get()->tag() << std::endl;
     t0_[le] = -1;
     i0_[le] = -1;
     t0_[le1] = -1;
     i0_[le1] = -1;        
   }
   
   void LinkMSD::update(const LinkRemoveEvent& event)
   {
     int le = (event.get()->tag())*2;
     int le1 = le+1;
     //std::cout << "RemoveEvent Tag: " << event.get()->tag() << std::endl;
     t0_[le] = -1;
     i0_[le] = -1;
     t0_[le1] = -1;
     i0_[le1] = -1;          
   }

   /*
   * Output param file after simulation is completed.
   */
   void LinkMSD::output() 
   {  
      // If outputFile_ was used to write block averages, close it.
      if (nSamplePerBlock_ != 0) {
         outputFile_.close();
      }
      
      // Write parameters to file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_); 
      outputFile_.close();
      
      // Write average to file
      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      int i;
      accumulator_[0].sample(0, outputFile_);;
      accumulator_[0].output(outputFile_);
      for(i=1; i<capacity_; ++i){
        accumulator_[i].output(outputFile_);
      } 
      outputFile_.close();
   }

}
#endif 
