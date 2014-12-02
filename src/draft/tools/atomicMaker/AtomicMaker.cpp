#include   "AtomicMaker.h"

#include <string>

void AtomicMaker::readParam(std::istream& in)
{
   read<Boundary>(in, "boundary", boundary_);
   readParamComposite(in, random_);
   read<int>(in, "nMolecule", nMolecule_);
}

void AtomicMaker::writeConfig(std::ostream& out)
{
   Vector r;
   Vector v; 
   int    iMol;

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
      boundary_.randomPosition(random_, r);
      out << r << std::endl;
      out << std::endl;
   }

}

int main() {
   AtomicMaker obj;
   obj.readParam(std::cin);
   obj.writeConfig(std::cout);
}
