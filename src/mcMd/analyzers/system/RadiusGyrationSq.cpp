/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RadiusGyrationSq.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/boundary/Boundary.h>
#include <simp/species/Species.h>
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   RadiusGyrationSq::RadiusGyrationSq(System& system) 
    : AverageAnalyzer<System>(system),
      speciesId_(-1)
   {  setClassName("RadiusGyrationSq"); }

   /*
   * Read parameters, allocate memory, and initialize accumulator.
   */
   void RadiusGyrationSq::readParameters(std::istream& in) 
   {
      AverageAnalyzer<System>::readParameters(in);
      read<int>(in, "speciesId", speciesId_);

      UTIL_CHECK(speciesId_ >= 0);
      UTIL_CHECK(speciesId_ < system().simulation().nSpecies());
      speciesPtr_ = &system().simulation().species(speciesId_);
      nAtom_ = speciesPtr_->nAtom();

      positions_.allocate(nAtom_); 
   }

   /*
   * Load state from an archive.
   */
   void RadiusGyrationSq::loadParameters(Serializable::IArchive& ar)
   {
      AverageAnalyzer::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      ar & nAtom_;
      ar & positions_;

      // Validate
      UTIL_CHECK(speciesId_ >= 0);
      UTIL_CHECK(speciesId_ < system().simulation().nSpecies());
      speciesPtr_ = &system().simulation().species(speciesId_);
      if (nAtom_ != speciesPtr_->nAtom()) {
         UTIL_THROW("Inconsistent values for nAtom");
      }
   }

   /*
   * Save state to archive.
   */
   void RadiusGyrationSq::save(Serializable::OArchive& ar)
   { ar & *this; }

   /* 
   * Evaluate squared radii of gyration of all chains, add to ensemble.
   */
   void RadiusGyrationSq::compute() 
   { 
      Molecule* moleculePtr;
      Vector    r1, r2, dR, Rcm;
      double    dRSq;
      int       i, j, nMolecule;

      dRSq = 0.0;
      nMolecule = system().nMolecule(speciesId_);
      for (i = 0; i < system().nMolecule(speciesId_); i++) {
         moleculePtr = &system().molecule(speciesId_, i);

         // Construct unwrapped map of molecule (no periodic b.c.'s)
         positions_[0] = moleculePtr->atom(0).position();
         Rcm = positions_[0];
         for (j = 1 ; j < nAtom_; j++) {
            r1 = moleculePtr->atom(j-1).position();
            r2 = moleculePtr->atom(j).position();
            system().boundary().distanceSq(r1, r2, dR);
            positions_[j]  = positions_[j-1];
            positions_[j] += dR;
            Rcm += positions_[j];
         }
         Rcm /= double(nAtom_);

         // Calculate dRSq
         for (j = 0 ; j < nAtom_; j++) {
            dR.subtract(positions_[j], Rcm);
            dRSq += dR.square();
         }
      }
      dRSq /= double(nMolecule);
      dRSq /= double(nAtom_);
      value_ = dRSq;
   }

}
