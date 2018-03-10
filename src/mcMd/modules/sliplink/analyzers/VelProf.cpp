#ifndef VEL_PROF_CPP
#define VEL_PROF_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VelProf.h"
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/boundary/Boundary.h>
#include <util/misc/FileMaster.h>        


namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /// Constructor.
   VelProf::VelProf(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      accumulator_(),
      totalMomentum_(0.0),
      velx_(),
      nSpec_(0),
      slabWidth_(0.0)
   {  setClassName("VelProf"); }

   /// Read parameters from file, and allocate data arrays.
   void VelProf::readParameters(std::istream& in) 
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"numberOfSlabs", nSlabs_);
      if (nSlabs_ <= 0) {
         UTIL_THROW("Invalid number of slabs");
      }
      
      //Allocate acumulator_ array.
      accumulator_. allocate(nSlabs_);
      
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      
      int i;
      for (i=0; i < nSlabs_; ++i){
        accumulator_[i].setNSamplePerBlock(nSamplePerBlock_);
      }
      // If nSamplePerBlock != 0, open an output file for block averages.
      if (nSamplePerBlock_ != 0) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }      
      
      // Allocate arrays 
      velx_.allocate(nSlabs_);
      nAtomi_.allocate(nSlabs_);     
   }

   
   /// Empty.
   void VelProf::setup() 
   {       
     nSpec_ = system().simulation().nSpecies(); 
     slabWidth_ = system().boundary().lengths()[2]/(double)nSlabs_;
     totalMomentum_ = 0.0;
     int i;
     for (i=0; i < nSlabs_; ++i){
        accumulator_[i].clear();
     }
     iBlock_ = 0;
   }
 
   /// Evaluate and output x-velocity profile as a function of z.
   void VelProf::sample(long iStep) 
   { 
     if (isAtInterval(iStep))  {
      int i, iSpec;
      System::MoleculeIterator  molIter;
      Molecule::AtomIterator  atomIter;
      double Lz = system().boundary().lengths()[2];
      double zcoord;
      int    imed=int(Lz/(2.0*slabWidth_));
      double vx, maxvx=-10.0, minvx=10.0, deltaP;
      Atom*  atomPtrmax=0;
      Atom*  atomPtrmin=0;
      
      for (i=0; i < nSlabs_; ++i){
        velx_[i] = 0.0;
	nAtomi_[i] = 0;
      }

      // Loop over all atoms
      for (iSpec = 0; iSpec < nSpec_; ++iSpec) { 
        for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter) {               
          for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {    
	    zcoord = (atomIter->position())[2]; 
	    i = int(zcoord/slabWidth_);	    
	    if(zcoord == Lz) i = nSlabs_-1;
	    vx = (atomIter->velocity())[0];
	    velx_[i] += vx;
	    nAtomi_[i] += 1;
	    if(iStep > 0){
	      if (i == 0) {
		if (vx < minvx){  
		  minvx = vx;
		  atomPtrmin = &(*atomIter);
		}
	      }
	      if (i == imed){
		if (vx > maxvx){  
		  maxvx = vx;
		  atomPtrmax = &(*atomIter);
		}
	      }
	    }  
	  } 
	}
      }
      
      // Calculate mean velocity for each slab.
      if(iStep > 0){
	iBlock_ += 1;
	
	for (i=0; i < nSlabs_; ++i){
	  if (nAtomi_[i] != 0 ){
	    velx_[i] = velx_[i]/(double)nAtomi_[i]; 
	  }
	  else {
	    velx_[i] = 0.0;
	  }
	  accumulator_[i].sample(velx_[i], outputFile_);
	}
	
	if (iBlock_ == nSamplePerBlock_) {
	  outputFile_ << std::endl;
	  iBlock_ = 0;
	}

        if(atomPtrmin!=0 && atomPtrmax!=0){
          (atomPtrmin->velocity())[0] = maxvx;
          (atomPtrmax->velocity())[0] = minvx;
	  deltaP = maxvx-minvx;
	  totalMomentum_ += deltaP;
	}
      } 
      
     } // if isAtInterval

   }

   /*
   * Output param file after simulation is completed.
   */
   void VelProf::output() 
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
      for(i=0; i<nSlabs_; ++i){
        accumulator_[i].output(outputFile_);
      }
      double A = system().boundary().lengths()[0] * system().boundary().lengths()[1];
      double nSteps = (double)system().simulation().iStep();
      outputFile_ << std::endl;
      outputFile_ << "Total transferred momentum: " << totalMomentum_ << std::endl;
      outputFile_ << "Total number of steps: " << nSteps << std::endl;
      outputFile_ << "Momentum flux: " << totalMomentum_/(2.0*A*nSteps) << std::endl;
      outputFile_.close();
   }

}
#endif 
