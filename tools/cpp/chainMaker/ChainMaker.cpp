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
   
   in >> outputStyle_;
}


void ChainMaker::writeChains(std::ostream& out)
{
   if (outputStyle_ == "McMd") {
      writeChainsMcMd(out);
   } else 
   if (outputStyle_ == "DdMd") {
      writeChainsDdMd(out);
   } else
   if (outputStyle_ == "DdMdMolecule") {
      writeChainsDdMdMole(out);
   } else {
      std::cout << "Unrecognized style " 
                << outputStyle_ << std::endl;
   }
}

void ChainMaker::writeChainsMcMd(std::ostream& out)
{
   Vector r;
   Vector v; 
   double beta = 1.0;
   int  bondType = 0;
   int  iMol, iAtom;

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
   int atomType = 0;
   int bondType = 0;
   int iMol, iAtom, i, j;

   out << "BOUNDARY" << std::endl;
   out << std::endl;
   out << boundary_ << std::endl;
   out << std::endl;
   out << "ATOMS" << std::endl;
   out << "nAtom  " << nMolecule_*nAtomPerMolecule_ << std::endl;

   i = 0;
   velocity.zero();
   for (iMol = 0; iMol < nMolecule_; ++iMol) {

      boundary_.randomPosition(random_, r);

      for (iAtom = 0; iAtom < nAtomPerMolecule_; ++iAtom) {

         out << Int(i,6) << Int(atomType, 10);
         for (j = 0; j < Dimension; ++j) {
            out << Dbl(r[j], 15, 6);
         }
         out << "  ";
         for (j = 0; j < Dimension; ++j) {
            out << Dbl(0.0, 12, 4);
         }
         out << std::endl;

         if (iAtom < nAtomPerMolecule_ - 1) {
            random_.unitVector(v);
            v *= bondPotential_.randomBondLength(&random_, beta, bondType);
            r += v;
            boundary_.shift(r);
         }

         ++i;
      }
   }

   // Write bonds
   out << std::endl;
   out << "BONDS" << std::endl;
   out << "nBond  " << nMolecule_*(nAtomPerMolecule_ -1 ) << std::endl;
   i = 0;
   j = 0;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {
      for (iAtom = 0; iAtom < nAtomPerMolecule_ - 1; ++iAtom) {
         out << Int(j,5) <<  Int(bondType, 5) << "  ";
         out << Int(i, 10) << Int(i + 1, 10) << std::endl;
         ++i;
         ++j;
      }
      ++i;
   }

   // Write angles
   int angleType = 0;
   out << std::endl;
   out << "ANGLES" << std::endl;
   out << "nAngle  " << nMolecule_*(nAtomPerMolecule_ -2) << std::endl;
   i = 0;
   j = 0;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {
      for (iAtom = 0; iAtom < nAtomPerMolecule_ - 2; ++iAtom) {
         out << Int(j,5) <<  Int(angleType, 5) << "  ";
         out << Int(i, 10) << Int(i + 1, 10) << Int(i + 2, 10)
             << std::endl;
         ++i;
         ++j;
      }
      i += 2;
   }

   // Write dihedrals
   int dihedralType = 0;
   out << std::endl;
   out << "DIHEDRALS" << std::endl;
   out << "nDihedral  " << nMolecule_*(nAtomPerMolecule_ -3) << std::endl;
   i = 0;
   j = 0;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {
      for (iAtom = 0; iAtom < nAtomPerMolecule_ - 3; ++iAtom) {
         out << Int(j,5) <<  Int(dihedralType, 5) << "  ";
         out << Int(i, 10) << Int(i + 1, 10) << Int(i + 2, 10)
             << Int(i + 3, 10) << std::endl;
         ++i;
         ++j;
      }
      i += 3;
   }

}

void ChainMaker::writeChainsDdMdMole(std::ostream& out)
{
   Vector r;
   Vector velocity;
   Vector v; 
   double beta = 1.0;
   int atomType = 0;
   int sId = 0;
   int bondType = 0;
   int iMol, iAtom, i, j;

   out << "BOUNDARY" << std::endl;
   out << std::endl;
   out << boundary_ << std::endl;
   out << std::endl;
   out << "ATOMS" << std::endl;
   out << "nAtom  " << nMolecule_*nAtomPerMolecule_ << std::endl;

   i = 0;
   velocity.zero();
   for (iMol = 0; iMol < nMolecule_; ++iMol) {

      boundary_.randomPosition(random_, r);

      for (iAtom = 0; iAtom < nAtomPerMolecule_; ++iAtom) {

         out << Int(i,10) << Int(atomType, 6)
             << Int(sId,6) << Int(iMol, 10) << Int(iAtom,6);
         out << "\n";
         for (j = 0; j < Dimension; ++j) {
            out << Dbl(r[j], 15, 6);
         }
         out << "\n";
         for (j = 0; j < Dimension; ++j) {
            out << Dbl(0.0, 15, 6);
         }
         out << "\n";

         if (iAtom < nAtomPerMolecule_ - 1) {
            random_.unitVector(v);
            v *= bondPotential_.randomBondLength(&random_, beta, bondType);
            r += v;
            boundary_.shift(r);
         }

         ++i;
      }
      out << std::endl;
   }

   // Write bonds
   out << std::endl;
   out << "BONDS" << std::endl;
   out << "nBond  " << nMolecule_*(nAtomPerMolecule_ -1 ) << std::endl;
   i = 0;
   j = 0;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {
      for (iAtom = 0; iAtom < nAtomPerMolecule_ - 1; ++iAtom) {
         out << Int(j,5) <<  Int(bondType, 5) << "  ";
         out << Int(i, 10) << Int(i + 1, 10) << std::endl;
         ++i;
         ++j;
      }
      ++i;
   }

   // Write angles
   int angleType = 0;
   out << std::endl;
   out << "ANGLES" << std::endl;
   out << "nAngle  " << nMolecule_*(nAtomPerMolecule_ -2) << std::endl;
   i = 0;
   j = 0;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {
      for (iAtom = 0; iAtom < nAtomPerMolecule_ - 2; ++iAtom) {
         out << Int(j,5) <<  Int(angleType, 5) << "  ";
         out << Int(i, 10) << Int(i + 1, 10) << Int(i + 2, 10)
             << std::endl;
         ++i;
         ++j;
      }
      i += 2;
   }

   // Write dihedrals
   int dihedralType = 0;
   out << std::endl;
   out << "DIHEDRALS" << std::endl;
   out << "nDihedral  " << nMolecule_*(nAtomPerMolecule_ -3) << std::endl;
   i = 0;
   j = 0;
   for (iMol = 0; iMol < nMolecule_; ++iMol) {
      for (iAtom = 0; iAtom < nAtomPerMolecule_ - 3; ++iAtom) {
         out << Int(j,5) <<  Int(dihedralType, 5) << "  ";
         out << Int(i, 10) << Int(i + 1, 10) << Int(i + 2, 10)
             << Int(i + 3, 10) << std::endl;
         ++i;
         ++j;
      }
      i += 3;
   }

}

int main() 
{
   ChainMaker obj;
   obj.readParam(std::cin);
   obj.writeChains(std::cout);
}

#endif
