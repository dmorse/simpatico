#include  <util/param/ParamComposite.h>
#include  <util/random/Random.h>
#include  <util/space/Vector.h>              // template parameter
#include  <util/containers/DArray.h>         // member template
#include  <util/boundary/Boundary.h>
#include  <simp/bond/HarmonicL0Bond.h>

#include <iostream>
#include <string>

using namespace Util;

class RingMaker : public ParamComposite
{

public:
   void readParam(std::istream& in);

   void writeRings(std::ostream& out);

private:

   Boundary       boundary_;
   Random         random_;
   Inter::HarmonicL0Bond bondPotential_;
   DArray<Vector> v_;
   int            nAtom_;
   int            nMolecule_;

};

void RingMaker::readParam(std::istream& in)
{
   //bondPotential_.setBoundary(boundary_);
   bondPotential_.setNBondType(1);

   read<Boundary>(in, "boundary", boundary_);
   readParamComposite(in, bondPotential_);
   readParamComposite(in, random_);
   read<int>(in, "nMolecule", nMolecule_);
   read<int>(in, "nAtomPerMolecule", nAtom_);

   if (nAtom_ < 2) UTIL_THROW("No. of atoms too small");
   v_.allocate(nAtom_);
}

void RingMaker::writeRings(std::ostream& out)
{
   Vector r, a, b;
   double beta = 1.0, length;
   int    bondType = 0;
   int    iMol, iAtom;

   out << "BOUNDARY" << std::endl;
   out << std::endl;
   out << "lengths " << boundary_ << std::endl;
   out << std::endl;
   out << "MOLECULES" << std::endl;
   out << std::endl;
   out << "species    " << 0 << std::endl;
   out << "nMolecule  " << nMolecule_ << std::endl;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {
      out << std::endl;
      out << "molecule   " << iMol << std::endl;

      // Generate the trial bond vectors
      v_[nAtom_-1].zero();
      for (iAtom = 1; iAtom < nAtom_; ++iAtom) {
         random_.unitVector(v_[iAtom-1]);
         v_[iAtom-1] *= 
           bondPotential_.randomBondLength(&random_, beta, bondType);
         v_[nAtom_-1] += v_[iAtom-1];
      }
      length = bondPotential_.randomBondLength(&random_, beta, bondType);
      length = v_[nAtom_-1].abs() - length;
      if (length > 0.0) {
         length /= v_[nAtom_-1].abs() * double(nAtom_-1);
         v_[nAtom_-1] *= length;
      } else {
         v_[nAtom_-1].zero();
      }

      // Choose first atom at random
      boundary_.randomPosition(random_, r);

      // Output the atom positions
      out << r << std::endl;
      for (iAtom = 1; iAtom < nAtom_; ++iAtom) {
         v_[iAtom-1] -= v_[nAtom_-1];
         r += v_[iAtom-1];
         boundary_.shift(r);
         out << r << std::endl;
      }

   }

}

// main
int main() {
   RingMaker ring;
   ring.readParam(std::cin);
   ring.writeRings(std::cout);
}
