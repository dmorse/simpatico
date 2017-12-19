#ifndef END_TO_END_XYZ_CPP
#define END_TO_END_XYZ_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EndtoEndXYZ.h"
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
   EndtoEndXYZ::EndtoEndXYZ(System& system) 
    : SystemAnalyzer<System>(system)
   {  setClassName("EndtoEndXYZ"); }

   /// Read parameters from file, and allocate data array.
   void EndtoEndXYZ::readParameters(std::istream& in) 
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulatorX_.setNSamplePerBlock(nSamplePerBlock_);
      accumulatorY_.setNSamplePerBlock(nSamplePerBlock_);
      accumulatorZ_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (nSamplePerBlock_ != 0) {
        fileMaster().openOutputFile(outputFileName("X.dat"), outputFileX_); 
	fileMaster().openOutputFile(outputFileName("Y.dat"), outputFileY_);
	fileMaster().openOutputFile(outputFileName("Z.dat"), outputFileZ_);
      }

      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }

      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      speciesPtr_ = &system().simulation().species(speciesId_);
      nAtom_ = speciesPtr_->nAtom();

      // Allocate an array of separation Vectors
      positions_.allocate(nAtom_); 
   }

   /*
   * Clear accumulator.
   */
   void EndtoEndXYZ::setup() 
   {  accumulatorX_.clear(); 
      accumulatorY_.clear();
      accumulatorZ_.clear();
   }
 
   /// Evaluate end-to-end vectors of all chains.
   void EndtoEndXYZ::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {

         Molecule* moleculePtr;
         Vector    r1, r2, dR;
         double    dxSq, dySq, dzSq;
         int       i, j, nMolecule;

         dxSq = 0.0;
         dySq = 0.0; 
         dzSq = 0.0; 
         nMolecule = system().nMolecule(speciesId_);
         for (i = 0; i < system().nMolecule(speciesId_); i++) {
            moleculePtr = &system().molecule(speciesId_, i);

            // Construct map of molecule with no periodic boundary conditions
            positions_[0] = moleculePtr->atom(0).position();
            for (j = 1 ; j < nAtom_; j++) {
               r1 = moleculePtr->atom(j-1).position();
               r2 = moleculePtr->atom(j).position();
               system().boundary().distanceSq(r1, r2, dR);
               positions_[j]  = positions_[j-1];
               positions_[j] += dR;
            }

            dR.subtract(positions_[0], positions_[nAtom_-1]);
            dxSq += dR[0]*dR[0];
            dySq += dR[1]*dR[1];
            dzSq += dR[2]*dR[2];    
            
     
         }
         dxSq /= double(nMolecule);
         dySq /= double(nMolecule);
         dzSq /= double(nMolecule); 
      
         accumulatorX_.sample(dxSq, outputFileX_);
         accumulatorY_.sample(dySq, outputFileY_); 
         accumulatorZ_.sample(dzSq, outputFileZ_); 

      } // if isAtInterval

   }

   /*
   * Output results to file after simulation is completed.
   */
   void EndtoEndXYZ::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (nSamplePerBlock_ != 0) {
         outputFileX_.close();
	 outputFileY_.close();
	 outputFileZ_.close();
      }

      // Write parameters to file
      fileMaster().openOutputFile(outputFileName("X.prm"), outputFileX_);
      ParamComposite::writeParam(outputFileX_); 
      outputFileX_.close();
      fileMaster().openOutputFile(outputFileName("Y.prm"), outputFileY_);
      ParamComposite::writeParam(outputFileY_); 
      outputFileY_.close();
      fileMaster().openOutputFile(outputFileName("Z.prm"), outputFileZ_);
      ParamComposite::writeParam(outputFileZ_); 
      outputFileZ_.close();      

      // Write average to file
      fileMaster().openOutputFile(outputFileName("X.ave"), outputFileX_);
      accumulatorX_.output(outputFileX_); 
      outputFileX_.close();
      fileMaster().openOutputFile(outputFileName("Y.ave"), outputFileY_);
      accumulatorY_.output(outputFileY_); 
      outputFileY_.close();
      fileMaster().openOutputFile(outputFileName("Z.ave"), outputFileZ_);
      accumulatorZ_.output(outputFileZ_); 
      outputFileZ_.close();      
   }

}
#endif 
