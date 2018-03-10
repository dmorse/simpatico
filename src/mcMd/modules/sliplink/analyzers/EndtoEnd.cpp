#ifndef END_TO_END_CPP
#define END_TO_END_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EndtoEnd.h"
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
   EndtoEnd::EndtoEnd(System& system) 
    : SystemAnalyzer<System>(system)
   {  setClassName("EndtoEnd"); }

   /// Read parameters from file, and allocate data array.
   void EndtoEnd::readParameters(std::istream& in) 
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
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
   void EndtoEnd::setup() 
   {  accumulator_.clear(); }
 
   /// Evaluate end-to-end vectors of all chains.
   void EndtoEnd::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {

         Molecule* moleculePtr;
         Vector    r1, r2, dR;
         double    dRSq;
         int       i, j, nMolecule;

         dRSq = 0.0;
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
            dRSq += dR.square();
            
     
         }
         dRSq /= double(nMolecule);
      
         accumulator_.sample(dRSq, outputFile_);

      } // if isAtInterval

   }

   /*
   * Output results to file after simulation is completed.
   */
   void EndtoEnd::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      // Write parameters to file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_); 
      outputFile_.close();

      // Write average to file
      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();
   }

}
#endif 
