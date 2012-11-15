#ifndef MCMD_RDF_CPP
#define MCMD_RDF_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "RDF.h"
#include <mcMd/simulation/Simulation.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/math/feq.h>
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>


#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   RDF::RDF(System& system) 
    : SystemDiagnostic<System>(system),
      outputFile_(),
      accumulator_(),
      typeNumbers_(),
      selector_(),
      max_(1.0),
      normSum_(0.0),
      nBin_(1),
      nAtomType_(0),
      isInitialized_(false)
   {  setClassName("RDF"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void RDF::readParameters(std::istream& in) 
   {
      Diagnostic::readParameters(in);
      read<double>(in, "max", max_);
      read<int>(in, "nBin", nBin_);
      read<PairSelector>(in, "selector", selector_);

      nAtomType_ = system().simulation().nAtomType();
      typeNumbers_.allocate(nAtomType_);
      accumulator_.setParam(max_, nBin_);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void RDF::loadParameters(Serializable::IArchive& ar)
   {
      Diagnostic::loadParameters(ar);
      loadParameter<double>(ar, "max", max_);
      loadParameter<int>(ar, "nBin", nBin_);
      loadParameter<PairSelector>(ar, "selector", selector_);

      ar & accumulator_;
      ar & nAtomType_;
      ar & typeNumbers_;
      ar & normSum_;

      // Validate
      if (nAtomType_ != system().simulation().nAtomType()) {
         UTIL_THROW("Inconsistent nAtomType");
      }
      if (nAtomType_ != typeNumbers_.capacity()) {
         UTIL_THROW("Inconsistent typeNumbers capacity");
      }
      if (nBin_ != accumulator_.nBin()) {
         UTIL_THROW("Inconsistent nBin values");
      }
      if (!feq(max_, accumulator_.max())) {
         UTIL_THROW("Inconsistent max values");
      }

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void RDF::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Add particle pairs to RDF histogram.
   */
   void RDF::setup() 
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

}
#endif
