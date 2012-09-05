#ifndef MCMD_BOND_LENGTH_DIST_CPP
#define MCMD_BOND_LENGTH_DIST_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BondLengthDist.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/util/FileMaster.h>
#include <util/archives/Serializable_includes.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   BondLengthDist::BondLengthDist(System& system) 
    : SystemDiagnostic<System>(system),
      isInitialized_(false)
   {  setClassName("BondLengthDist"); }

   /// Read parameters from file, and allocate data array.
   void BondLengthDist::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      readParamComposite(in, accumulator_);
      accumulator_.clear();
      isInitialized_ = true;
   }

   /*
   * Clear accumulator.
   */
   void BondLengthDist::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Object is not initialized"); 
      }
      accumulator_.clear();
   }
 
   /// Add particle pairs to RDF histogram.
   void BondLengthDist::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         System::MoleculeIterator molIter;
         Molecule::BondIterator   bondIter;
         double    lsq, l;

         for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(bondIter); bondIter.notEnd(); ++bondIter) {
               lsq = system().boundary().distanceSq( bondIter->atom(0).position(), 
                                                     bondIter->atom(1).position());
               l = sqrt(lsq);
               accumulator_.sample(l);
            }
         }
      }
   }  


   /// Output results to file after simulation is completed.
   void BondLengthDist::output() 
   {  

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

   /*
   * Save state to binary file archive.
   */
   void BondLengthDist::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void BondLengthDist::load(Serializable::IArchiveType& ar)
   { ar & *this; }
}
#endif
