#ifndef G1MSD_CPP
#define G1MSD_CPP

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "G1MSD.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>
#include <simp/boundary/Boundary.h>
#include <util/space/Dimension.h>

namespace McMd
{
  
   using namespace Util; 
   using namespace Simp; 
  
   /*
   * Constructor.
   */
   G1MSD::G1MSD(System& system) :
      SystemAnalyzer<System>(system),
      outputFile_(),
      accumulator_(),
      truePositions_(),
      oldPositions_(),
      shifts_(),
      speciesPtr_(0),
      speciesId_(-1),
      nMolecule_(-1),
      capacity_(-1),
      nAtom_(0)
   {  setClassName("G1MSD"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void G1MSD::readParameters(std::istream& in) 
   {

      // Read interval and parameters for AutoCorrArray
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "capacity", capacity_);

      // Validate parameters
      if (speciesId_ < 0)       UTIL_THROW("Negative speciesId");
      if (speciesId_ >= system().simulation().nSpecies()) 
                                UTIL_THROW("speciesId > nSpecies");
      if (capacity_ <= 0)       UTIL_THROW("Negative capacity");

      speciesPtr_ = &system().simulation().species(speciesId_);
      nAtom_ = speciesPtr_->nAtom();
   }

   /*
   * Allocate arrays.
   */
   void G1MSD::setup() 
   {
 
      // Get number of molecules of this species in the System. 
      nMolecule_ = system().nMolecule(speciesId_);
      
      int nratoms;
      nratoms =  nMolecule_*nAtom_;
      // Allocate arrays of position Vectors and shifts
      truePositions_.allocate(nratoms); 
      oldPositions_.allocate(nratoms); 
      shifts_.allocate(nratoms); 

      // Initialize the AutoCorrArray object
      accumulator_.setParam(nratoms, capacity_);

      // Store initial positions, and set initial shift vectors.
      Vector     r;
      IntVector  zero(0);
      Molecule*  moleculePtr;
      int        iatom = 0;
      for (int i = 0; i < nMolecule_; ++i) {
	 moleculePtr = &system().molecule(speciesId_, i);
	 for (int j = 0 ; j < nAtom_; j++) {
            r = moleculePtr->atom(j).position();
            system().boundary().shift(r);
            oldPositions_[iatom] = r;
            shifts_[iatom] = zero;
	    iatom++;
	 }
      }

   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void G1MSD::sample(long iStep) 
   { 
      if (!isAtInterval(iStep)) return;

      // Confirm that nMolecule has remained constant, and nMolecule > 0.
      if (nMolecule_ <= 0) {
         UTIL_THROW("nMolecule <= 0");
      }
      if (nMolecule_ != system().nMolecule(speciesId_)) {
         UTIL_THROW("Number of molecules has changed.");
      }

      Vector     r;
      IntVector  shift;
      Molecule*  moleculePtr;
      Vector     lengths = system().boundary().lengths();
      int        i, j, k;
      int        iatom = 0;
      for (i = 0; i < nMolecule_; ++i) {
	 moleculePtr = &system().molecule(speciesId_, i);
	 for (j = 0 ; j < nAtom_; j++) {	 
            r = moleculePtr->atom(j).position();
            system().boundary().shift(r);

            // Compare current r to previous position, oldPositions_[iatom]
            system().boundary().distanceSq(r, oldPositions_[iatom], shift);

            // If this atom crossed a boundary, increment its shift vector
            shifts_[iatom] += shift;

            // Reconstruct true position
            for (k = 0; k < Dimension; ++k) {
               truePositions_[iatom][k] = r[k] + shifts_[iatom][k]*lengths[k];
            }

            // Store current position in box for comparison to next one
            oldPositions_[iatom] = r;
	    
	    iatom++;
	 }
      }
      accumulator_.sample(truePositions_);

   }

   /// Output results to file after simulation is completed.
   void G1MSD::output() 
   {  

      // Echo parameter to analyzer log file
      fileMaster().openOutputFile(outputFileName(), outputFile_);
      writeParam(outputFile_); 
      outputFile_ << std::endl;
      outputFile_ << "nrAtoms  " << accumulator_.nEnsemble() << std::endl;
      outputFile_ << "buffercapacity " << accumulator_.bufferCapacity()  << std::endl;
      outputFile_ << "nSample    " << accumulator_.nSample()   << std::endl;
      outputFile_ << std::endl;
      outputFile_ << "Format of *.dat file:" << std::endl;
      outputFile_ << "[time in samples]  [Mean Sq Displacement]"
                  << std::endl;
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}

#endif 
