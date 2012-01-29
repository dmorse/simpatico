#ifndef CHAIN_MAKER_CPP
#define CHAIN_MAKER_CPP

#include   "ChainMaker.h"

#include <string>

void ChainMaker::readParam(std::istream& in)
{
   bondPotential_.setNBondType(1);

   read<Boundary>(in, "boundary", boundary_);
   readParamComposite(in, bondPotential_);
   readParamComposite(in, random_);
   read<int>(in, "nMolecule", nMolecule_);
   read<int>(in, "nAtomPerMolecule", nAtomPerMolecule_);
}

void ChainMaker::writeChains(std::ostream& out)
{
   Vector r;
   Vector v; 
   double beta = 1.0;
   int    bondType = 0;
   int    iMol, iAtom;

   out << "BOUNDARY" << std::endl;
   out << std::endl;
   out << boundary_ << std::endl;
   out << std::endl;
   out << "MOLECULES" << std::endl;
   out << std::endl;
   out << "species    " << 0 << std::endl;
   out << "nMolecule  " << nMolecule_ << std::endl;
   out << std::endl;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {
      out << "molecule   " << iMol << std::endl;
      // Choosen first atom at random
      boundary_.randomPosition(random_, r);
      out << r << std::endl;
      for (iAtom = 1; iAtom < nAtomPerMolecule_; ++iAtom) {
         random_.unitVector(v);
         v *= bondPotential_.randomBondLength(&random_, beta, bondType);
         r += v;
         boundary_.shift(r);
         out << r << std::endl;
      }
      out << std::endl;
   }

}

int main() {
   ChainMaker obj;
   obj.readParam(std::cin);
   obj.writeChains(std::cout);
}

#endif
