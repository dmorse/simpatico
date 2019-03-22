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
#include <util/archives/Serializable_includes.h>

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
      DistributionAnalyzer<System>::readParameters(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "typeId", typeId_);
      UTIL_CHECK(speciesId_ >= -1);
      UTIL_CHECK(typeId_ >= -1);
   }

   /*
   * Load state from an archive.
   */
   void DihedralDist::loadParameters(Serializable::IArchive& ar)
   {
      DistributionAnalyzer<System>::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "typeId", typeId_);

      // Validate
      UTIL_CHECK(speciesId_ < -1);
      UTIL_CHECK(typeId_ < -1);
   }

   /// Add particle pairs to RDF histogram.
   void DihedralDist::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         Boundary& boundary = system().boundary();
         Torsion torsion;
         System::MoleculeIterator mIter;
         Molecule::DihedralIterator dIter;  
         Vector dr1, dr2, dr3;
         double phi;

         system().begin(speciesId_, mIter);
         for ( ; mIter.notEnd(); ++mIter) {
            for ( mIter->begin(dIter); dIter.notEnd(); ++dIter) {
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
