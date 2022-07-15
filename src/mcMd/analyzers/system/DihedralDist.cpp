/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DihedralDist.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>

#include <simp/species/Species.h>
#include <simp/boundary/Boundary.h>
#include <simp/interaction/dihedral/Torsion.h>

#include <util/space/Vector.h>
#include <util/math/feq.h>
#include <util/misc/FileMaster.h>
//#include <util/archives/Serializable_includes.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /// Constructor.
   DihedralDist::DihedralDist(System& system) 
    : DistributionAnalyzer<System>(system),
      speciesId_(-1),
      typeId_(-1)
   {  setClassName("DihedralDist"); }

   /// Read parameters from file, and allocate data array.
   void DihedralDist::readParameters(std::istream& in) 
   {
      Analyzer::readParameters(in);
      min_ = 0.0;
      max_ = 2.0*acos(0.0);
      read<int>(in, "nBin", nBin_);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "typeId", typeId_);

      UTIL_CHECK(nBin_ > 0);
      UTIL_CHECK(speciesId_ >= -1);
      UTIL_CHECK(typeId_ >= -1);

      accumulator().setParam(min_, max_, nBin_);
      accumulator().clear();
   }

   /*
   * Load state from an archive.
   */
   void DihedralDist::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      ar >> min_;
      ar >> max_;
      loadParameter<int>(ar, "nBin", nBin_);
      ar >> accumulator();
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "typeId", typeId_);

      // Validate
      UTIL_CHECK(feq(min_, 0.0));
      UTIL_CHECK(feq(max_, 2.0*acos(0.0)));
      UTIL_CHECK(nBin_ > 0);
      UTIL_CHECK(feq(accumulator().min(), 0.0));
      UTIL_CHECK(feq(accumulator().min(), min_));
      UTIL_CHECK(feq(accumulator().max(), max_));
      UTIL_CHECK(accumulator().nBin() == nBin_);
      UTIL_CHECK(speciesId_ >= -1);
      UTIL_CHECK(typeId_ >= -1);

   }

   /*
   * Save state to an archive (implicitly calls serialize).
   */
   void DihedralDist::save(Serializable::OArchive& ar)
   {  ar & *this; }

   void DihedralDist::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         int nSpecies = system().simulation().nSpecies();
         UTIL_CHECK(speciesId_ >=-1);
         if (speciesId_ >= 0) {
            UTIL_CHECK(speciesId_ < nSpecies);
            sampleSpecies(speciesId_);
         } else {
            for (int is = 0; is < nSpecies; ++is) {
               sampleSpecies(is);
            }
         } 

      }
   }

   /*
   * Sample dihedrals for species with index is
   */
   void DihedralDist::sampleSpecies(int is) 
   {
      Boundary& boundary = system().boundary();
      Torsion torsion;
      System::MoleculeIterator mIter;
      Molecule::DihedralIterator dIter;  
      Vector dr1, dr2, dr3;
      double phi;

      system().begin(is, mIter);
      for ( ; mIter.notEnd(); ++mIter) {
         for ( mIter->begin(dIter); dIter.notEnd(); ++dIter) {
            if (typeId_ == -1 || typeId_ == dIter->typeId()) {
               boundary.distanceSq(dIter->atom(1).position(),
                                   dIter->atom(0).position(), dr1);
               boundary.distanceSq(dIter->atom(2).position(),
                                   dIter->atom(1).position(), dr2);
               boundary.distanceSq(dIter->atom(3).position(),
                                   dIter->atom(2).position(), dr3);
               torsion.computeAngle(dr1, dr2, dr3);
               phi = torsion.phi();
               increment(phi);
            }
         }
      }
   }

}
