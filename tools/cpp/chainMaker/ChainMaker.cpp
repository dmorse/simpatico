#ifndef CHAIN_MAKER_CPP
#define CHAIN_MAKER_CPP

#include   "ChainMaker.h"
#include   <util/format/Int.h>
#include   <util/format/Dbl.h>

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



void ChainMaker::writeChainsMcMd(std::ostream& out)
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

void ChainMaker::writeChainsDdMd(std::ostream& out)
{
   Vector r;
   Vector velocity;
   Vector v; 
   double beta = 1.0;
   int    bondType = 0;
   int    iMol, iAtom, i, j;

   out << "BOUNDARY" << std::endl;
   out << std::endl;
   out << boundary_ << std::endl;
   out << std::endl;
   out << "ATOMS" << std::endl;
   out << "nAtom  " << nMolecule_*nAtomPerMolecule_ << std::endl;
   out << std::endl;

   i = 0;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {

      boundary_.randomPosition(random_, r);
      ++i;

      out << Int(i,6) << Int(bondType, 5);
      out << r << std::endl;

      for (iAtom = 1; iAtom < nAtomPerMolecule_; ++iAtom) {
         random_.unitVector(v);
         v *= bondPotential_.randomBondLength(&random_, beta, bondType);
         r += v;
         boundary_.shift(r);

         out << Int(i,6) << Int(bondType, 5);
         out << r << std::endl;


         ++i;
      }
   }

   out << std::endl;
   out << "BONDS" << std::endl;
   out << "nBond  " << nMolecule_*(nAtomPerMolecule_ -1 ) << std::endl;
   out << std::endl;
   i = 0;
   j = 0;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {

      for (iAtom = 1; iAtom < nAtomPerMolecule_; ++iAtom) {
         out << Int(j,5) <<  Int(bondType, 5) << "  ";
         out << Int(i, 5) << Int(i + 1, 5) << std::endl;
         ++i;
         ++j;
      }
      ++i;
    }
}

int main() {
   ChainMaker obj;
   obj.readParam(std::cin);
   int format;
   std::cin >> format;
   if (format == 0) {
      obj.writeChainsMcMd(std::cout);
   } else 
   if (format == 1) {
      obj.writeChainsDdMd(std::cout);
   }
}

#endif
