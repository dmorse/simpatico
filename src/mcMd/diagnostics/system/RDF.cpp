#ifndef RDF_CPP
#define MCMD_RDF_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "RDF.h"
#include <mcMd/simulation/Simulation.h>
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
   RDF::RDF(System& system) 
    : SystemDiagnostic<System>(system),
      isInitialized_(false)
   {}

   /// Read parameters from file, and allocate data array.
   void RDF::readParam(std::istream& in) 
   {

      // Read interval and parameters for AutoCorrArray
      readInterval(in);
      readOutputFileName(in);
      readParamComposite(in, accumulator_);

      read<PairSelector>(in, "selector", selector_);

      nAtomType_ = system().simulation().nAtomType();
      typeNumbers_.allocate(nAtomType_);

      isInitialized_ = true;
   }

   /// Add particle pairs to RDF histogram.
   void RDF::initialize() 
   { 
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized");
      }
      accumulator_.clear(); 
   }

   /// Add particle pairs to RDF histogram.
   void RDF::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         accumulator_.beginSnapshot();

         System::ConstMoleculeIterator molIter1, molIter2;
         Molecule::ConstAtomIterator   atomIter1, atomIter2;
         Vector    r1, r2;
         double    dRsq, dR;
         Boundary* boundaryPtr;
         int       iSpecies1, iSpecies2, nSpecies, i;
   
         for (i = 0; i < nAtomType_; ++i) {
            typeNumbers_[i] = 0;
         }

         boundaryPtr = &system().boundary();
         nSpecies    = system().simulation().nSpecies();

         // Loop over atom 1
         for (iSpecies1 = 0; iSpecies1 < nSpecies; ++iSpecies1) {
            system().begin(iSpecies1, molIter1); 
            for ( ; molIter1.notEnd(); ++molIter1) {
               molIter1->begin(atomIter1); 
               for ( ; atomIter1.notEnd(); ++atomIter1) {
                  r1 = atomIter1->position();
 
                  ++typeNumbers_[atomIter1->typeId()];
  
                  // Loop over atom 2 
                  for (iSpecies2 = 0; iSpecies2 < nSpecies; ++iSpecies2) {
                     system().begin(iSpecies2, molIter2); 
                     for ( ; molIter2.notEnd(); ++molIter2) {

                        //Check if molecules are the same  
                        //if ( &(*molIter2) != &(*molIter1)) {

                           molIter2->begin(atomIter2);
                           for ( ; atomIter2.notEnd(); ++atomIter2) {

                              if (selector_.match(*atomIter1, *atomIter2)) {
                                 r2 = atomIter2->position();
   
                                 dRsq = boundaryPtr->distanceSq(r1, r2);
                                 dR   = sqrt(dRsq);

                                 accumulator_.sample(dR);

                              }
               
                           }

                        //}

                     }
                  } // for iSpecies2
   
               }
            }
         } // for iSpecies1

         // Increment normSum_
         double number = 0;
         for (i = 0; i < nAtomType_; ++i) {
            number  += typeNumbers_[i];
         }
         normSum_ += number*number/boundaryPtr->volume();

      } // if isAtInterval

   }


   /// Output results to file after simulation is completed.
   void RDF::output() 
   {  

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      double nSnapshot = double(accumulator_.nSnapshot()) ;
      accumulator_.setNorm(normSum_ / nSnapshot);
      if (selector_.pairType() == PairSelector::ALL) {
         accumulator_.setOutputIntegral(true);
      }
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

   /*
   * Save state to binary file archive.
   */
   void RDF::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void RDF::load(Serializable::IArchiveType& ar)
   { ar & *this; }
}
#endif
