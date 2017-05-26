/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsConfigIo.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>

#include <simp/species/Species.h>

#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/param/Label.h>
#include <util/format/Format.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>

#include <vector>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /* 
   * Constructor.   
   */
   LammpsConfigIo::LammpsConfigIo(System &system) 
   : ConfigIo(system)
   {}
 
   /* 
   * Destructor.   
   */
   LammpsConfigIo::~LammpsConfigIo() 
   {}

   void LammpsConfigIo::read(std::istream &in)
   {
      // Calculate atomCapacity for entire simulation
      int atomCapacity = 0;
      int bondCapacity = 0;
      int nSpecies = simulation().nSpecies();
      int speciesCapacity = 0;
      int iSpec, i;
      Species* speciesPtr;
      for (iSpec = 0; iSpec < nSpecies; ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         speciesCapacity = speciesPtr->capacity();
         atomCapacity += speciesCapacity*speciesPtr->nAtom();
         #ifdef SIMP_BOND
         bondCapacity += speciesCapacity*speciesPtr->nBond();
         #endif
      }


      // Read and discard title line
      std::string       line;
      std::getline(in, line);

      // Read numbers of atoms, bonds, etc.
      int nAtom;
      int nBond;
      int nAngle;
      int nDihedral;
      int nImproper;
      in >> nAtom >> Label("atoms");
      in >> nBond >> Label("bonds");
      in >> nAngle >> Label("angles");
      in >> nDihedral >> Label("dihedrals");
      in >> nImproper >> Label("impropers");
  
      /* 
      * Validate nAtom and nBond
      * Lammps files can be read only if the number of atoms and bonds
      * in the lammps file exactly matches the corresponding capacities.
      */
      if (nAtom != atomCapacity) {
         UTIL_THROW("nAtom != atomCapacity");
      } 
      if (nBond != bondCapacity) {
         UTIL_THROW("nAtom != atomCapacity");
      } 

      // Read numbers of atom types, bond types, etc.
      int nAtomType;
      int nBondType;
      int nAngleType;
      int nDihedralType;
      int nImproperType;
      in >> nAtomType >> Label("atom") >> Label("types");
      in >> nBondType >> Label("bond") >> Label("types");
      in >> nAngleType >> Label("angle") >> Label("types");
      in >> nDihedralType >> Label("dihedral") >> Label("types");
      in >> nImproperType >> Label("improper") >> Label("types");

      if (nAtomType > simulation().nAtomType()) {
         UTIL_THROW("nAtomType > simulation().nAtomType()");
      } 
      // Read boundary dimensions
      Vector lengths;
      Vector min;
      Vector max;
      in >> min[0] >> max[0] >> Label("xlo") >> Label("xhi");
      in >> min[1] >> max[1] >> Label("ylo") >> Label("yhi");
      in >> min[2] >> max[2] >> Label("zlo") >> Label("zhi");
      lengths.subtract(max, min);
      boundary().setOrthorhombic(lengths);

      // Read particle masses (discard values)
      double mass;
      int atomTypeId;
      in >> Label("Masses");
      for (i = 0; i < nAtomType; ++i) {
         in >> atomTypeId >> mass;
      }

      /*
      * Read atomic positions
      *
      * Atom tags must appear in order, numbered from 1
      */
      in >> Label("Atoms");
      Molecule* molPtr;
      Molecule::AtomIterator atomIter;
      IntVector shift;
      int atomId, molId, iMol, nMolecule;
      for (iSpec=0; iSpec < nSpecies; ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         nMolecule = speciesPtr->capacity();
         for (iMol = 0; iMol < nMolecule; ++iMol) {
            molPtr = &(simulation().getMolecule(iSpec));
            system().addMolecule(*molPtr);
   
            if (molPtr != &system().molecule(iSpec, iMol)) {
               UTIL_THROW("Molecule index error");
            }
   
            // Read positions, shift into primary cell
            for (molPtr->begin(atomIter); atomIter.notEnd(); ++atomIter) {

               in >> atomId >> molId >> atomTypeId;
               if (atomId != atomIter->id() + 1) 
                  UTIL_THROW("Atom tags not ordered");

               in >> atomIter->position();
 
               // Shift coordinates to box (0, length(i));
               atomIter->position() += min;

               in >> shift;

            }
         }
      }

   }

   void LammpsConfigIo::write(std::ostream &out)
   {

      // Write first line (skipped) and a blank line.
      out << "LAMMPS data file" << "\n";
      out << "\n";

      // Count total numbers of atoms and bonds in all species.
      Species* speciesPtr;
      int iSpec, nMolecule;
      int nAtom = 0;
      int nBond = 0;
      int nAngle = 0;
      int nDihedral = 0;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         nMolecule = system().nMolecule(iSpec);
         nAtom += nMolecule*(speciesPtr->nAtom());
         #ifdef SIMP_BOND
         nBond += nMolecule*(speciesPtr->nBond());
         #endif
         #ifdef SIMP_ANGLE
         nAngle += nMolecule*(speciesPtr->nAngle());
         #endif
         #ifdef SIMP_DIHEDRAL
         nDihedral += nMolecule*(speciesPtr->nDihedral());
         #endif
      }

      // Write numbers of atoms, bonds, etc.
      Format::setDefaultWidth(8);
      out << Int(nAtom)  << " atoms    " << "\n";
      out << Int(nBond)  << " bonds    " << "\n";
      out << Int(nAngle) << " angles   " << "\n";
      out << Int(nDihedral) << " dihedrals" << "\n";
      out << Int(0) << " impropers" << "\n";
      out << "\n";

      // Write numbers of atom types, bond types, etc.
      Format::setDefaultWidth(5);
      out << Int(simulation().nAtomType()) << " atom types" << "\n";
      int nBondType = 0;
      int nAngleType = 0;
      int nDihedralType = 0;
      #ifdef SIMP_BOND
      nBondType = simulation().nBondType();
      #endif
      #ifdef SIMP_ANGLE
      nAngleType = simulation().nAngleType();
      #endif
      #ifdef SIMP_DIHEDRAL
      nDihedralType = simulation().nDihedralType();
      #endif
      out << Int(nBondType) << " bond types" << "\n";
      out << Int(nAngleType) << " angle types" << "\n";
      out << Int(nDihedralType) << " dihedral types" << "\n";
      out << Int(0) << " improper types" << "\n";
      out << "\n";

      // Write Boundary dimensions
      Vector lengths = boundary().lengths();
      Format::setDefaultWidth(15);
      Format::setDefaultPrecision(7);
      out << Dbl(0.0) << Dbl(lengths[0]) << "  xlo xhi" << "\n";
      out << Dbl(0.0) << Dbl(lengths[1]) << "  ylo yhi" << "\n";
      out << Dbl(0.0) << Dbl(lengths[2]) << "  zlo zhi" << "\n";
      out << "\n";

      // Write masses (all set to 1.0 for now)
      // lammps atom type = Simpatico atom type + 1
      out << "Masses" << "\n";
      out << "\n";
      for (int iType = 0; iType < simulation().nAtomType(); ++iType) {
          out << Int(iType+1, 5) << Dbl(1.0) << "\n";
      }
      out << "\n";

      // Write atomic positions
      // lammps atom tag = Simpatico atom id + 1
      // lammps molecule id = Simpatico molecule id + 1
      out << "Atoms" << "\n";
      out << "\n";
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      int i;
      int shift = 0;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               out << Int( atomIter->id() + 1, 8 ) << Int( molIter->id() + 1, 8 );
               out << Int( atomIter->typeId() + 1, 5 );
               out << atomIter->position();
               for (i = 0; i < Dimension; ++i) {
                  out << Int(shift, 4);
               }
               out << "\n";
            }
         }
      }
      out << "\n";

      #ifdef SIMP_BOND
      if (nBond > 0) {
         out << "Bonds" << "\n";
         out << "\n";
         Molecule::BondIterator bondIter;
         int iBond = 1;
         for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
            for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               for (molIter->begin(bondIter); bondIter.notEnd(); ++bondIter) {
                  out << Int(iBond, 8 ) << Int(bondIter->typeId() + 1, 5);
                  out << Int(bondIter->atom(0).id() + 1, 8);
                  out << Int(bondIter->atom(1).id() + 1, 8);
                  out << "\n";
                  ++iBond;
               }
            }
         }
         out << "\n";
      }
      #endif

      #ifdef SIMP_ANGLE
      if (nAngle > 0) {
         out << "Angles" << "\n";
         out << "\n";
         Molecule::AngleIterator angleIter;
         int iAngle = 1;
         for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
            for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               for (molIter->begin(angleIter); angleIter.notEnd(); ++angleIter) {
                  out << Int(iAngle, 8) << Int(angleIter->typeId() + 1, 5);
                  out << Int(angleIter->atom(0).id() + 1, 8);
                  out << Int(angleIter->atom(1).id() + 1, 8);
                  out << Int(angleIter->atom(2).id() + 1, 8);
                  out << "\n";
                  ++iAngle;
               }
            }
         }
         out << "\n";
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      if (nDihedral > 0) {
         out << "Dihedrals" << "\n";
         out << "\n";
         Molecule::DihedralIterator dihedralIter;
         int iDihedral = 1;
         for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
            for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               molIter->begin(dihedralIter); 
               for ( ; dihedralIter.notEnd(); ++dihedralIter) {
                  out << Int(iDihedral, 8); 
                  out << Int(dihedralIter->typeId() + 1, 5);
                  out << Int(dihedralIter->atom(0).id() + 1, 8);
                  out << Int(dihedralIter->atom(1).id() + 1, 8);
                  out << Int(dihedralIter->atom(2).id() + 1, 8);
                  out << Int(dihedralIter->atom(3).id() + 1, 8);
                  out << "\n";
                  ++iDihedral;
               }
            }
         }
         out << "\n";
      }
      #endif

      // Reset Format defaults to initialization values
      Format::initStatic();
   }

} 
