#ifndef MCMD_DDMD_CONFIG_IO_CPP
#define MCMD_DDMD_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdConfigIo.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
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

   /* 
   * Constructor.   
   */
   DdMdConfigIo::DdMdConfigIo(System &system) 
   : ConfigIo(system)
   {}
 
   /* 
   * Destructor.   
   */
   DdMdConfigIo::~DdMdConfigIo() 
   {}

   void DdMdConfigIo::read(std::istream &in)
   {
      // Calculate atom and bond capacities for entire simulation
      int atomCapacity = 0;
      int bondCapacity = 0;
      int nSpecies = simulation().nSpecies();
      int speciesCapacity = 0;
      int iSpec;
      Species* speciesPtr;
      for (iSpec = 0; iSpec < nSpecies; ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         speciesCapacity = speciesPtr->capacity();
         atomCapacity   += speciesCapacity*speciesPtr->nAtom();
         bondCapacity   += speciesCapacity*speciesPtr->nBond();
      }

      // Read boundary
      in >> Label("BOUNDARY");
      in >> system().boundary();

      // Read atomic positions
      // Atom tags must appear in order, numbered from 0
      // Number of atoms must match atomCapacity.
      in >> Label("ATOMS");
      int nAtom;
      in >> Label("nAtom") >> nAtom;
      if (nAtom != atomCapacity) {
         UTIL_THROW("nAtom != atomCapacity");
      }
      Molecule*                molPtr;
      Molecule::AtomIterator   atomIter;
      int iMol, nMolecule, atomId, atomTypeId;
      for (iSpec=0; iSpec < nSpecies; ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         nMolecule  = speciesPtr->capacity();
         for (iMol = 0; iMol < nMolecule; ++iMol) {
            molPtr = &(speciesPtr->reservoir().pop());
            system().addMolecule(*molPtr);
   
            // Read positions.
            for (molPtr->begin(atomIter); atomIter.notEnd(); ++atomIter) {

               in >> atomId >> atomTypeId;
               if (atomId != atomIter->id()) {
                  UTIL_THROW("Atom tags not ordered");
               }

               in >> atomIter->position();
               in >> atomIter->velocity();
 
            }
         }
      }

   }

   void DdMdConfigIo::write(std::ostream &out)
   {
      // Count total numbers of atoms and bonds in all species.
      Species  *speciesPtr;
      int iSpec, nMolecule;
      int nAtom = 0;
      int nBond = 0;
      #ifdef INTER_ANGLE
      int nAngle = 0;
      #endif
      #ifdef INTER_DIHEDRAL
      int nDihedral = 0;
      #endif
      for (iSpec = 0; iSpec < simulation().nSpecies(); ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         nMolecule  = system().nMolecule(iSpec);
         nAtom += nMolecule*(speciesPtr->nAtom());
         nBond += nMolecule*(speciesPtr->nBond());
         #ifdef INTER_ANGLE
         nAngle += nMolecule*(speciesPtr->nAngle());
         #endif
         #ifdef INTER_DIHEDRAL
         nDihedral += nMolecule*(speciesPtr->nDihedral());
         #endif
      }

      // Write boundary
      out << "BOUNDARY" << std::endl;
      out << system().boundary() << std::endl;
      out << std::endl;

      // Write atoms
      out << "ATOMS" << std::endl;
      out << "nAtom  " << nAtom << std::endl;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      int i = 0;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         system().begin(iSpec, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            molIter->begin(atomIter); 
            for ( ; atomIter.notEnd(); ++atomIter) {
               out << Int(atomIter->id(), 10);
               out << Int(atomIter->typeId(), 5);
               out << atomIter->position();
               out << atomIter->velocity();
               out << std::endl;
               ++i;
            }
         }
      }
      out << std::endl;

      // Write Bonds
      out << "BONDS" << std::endl;
      out << "nBond  " << nBond << std::endl;
      Molecule::BondIterator bondIter;
      i = 0;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (system().simulation().species(iSpec).nBond() > 0) {
            system().begin(iSpec, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(bondIter); 
               for ( ; bondIter.notEnd(); ++bondIter) {
                  out << Int(i, 8) << Int(bondIter->typeId(), 5);
                  out << Int(bondIter->atom(0).id(), 10);
                  out << Int(bondIter->atom(1).id(), 10);
                  out << std::endl;
                  ++i;
               }
            }
         }
      }
      out << std::endl;

      #ifdef INTER_ANGLE
      // Write Angles
      out << "ANGLES" << std::endl;
      out << "nAngle  " << nAngle << std::endl;
      Molecule::AngleIterator angleIter;
      i = 0;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (system().simulation().species(iSpec).nAngle() > 0) {
            system().begin(iSpec, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(angleIter); 
               for ( ; angleIter.notEnd(); ++angleIter) {
                  out << Int(i, 8) << Int(angleIter->typeId(), 5);
                  out << Int(angleIter->atom(0).id(), 10);
                  out << Int(angleIter->atom(1).id(), 10);
                  out << Int(angleIter->atom(2).id(), 10);
                  out << std::endl;
                  ++i;
               }
            }
         }
      }
      out << std::endl;
      #endif

      #ifdef INTER_DIHEDRAL
      // Write Dihedral
      out << "DIHEDRALS" << std::endl;
      out << "nDihedral  " << nDihedral << std::endl;
      Molecule::DihedralIterator dihedralIter;
      i = 0;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (system().simulation().species(iSpec).nDihedral() > 0) {
            system().begin(iSpec, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(dihedralIter); 
               for ( ; dihedralIter.notEnd(); ++dihedralIter) {
                  out << Int(i, 8) << Int(dihedralIter->typeId(), 5);
                  out << Int(dihedralIter->atom(0).id(), 10);
                  out << Int(dihedralIter->atom(1).id(), 10);
                  out << Int(dihedralIter->atom(2).id(), 10);
                  out << Int(dihedralIter->atom(3).id(), 10);
                  out << std::endl;
                  ++i;
               }
            }
         }
      }
      out << std::endl;
      #endif

      // Reset Format defaults to initialization values
      Format::initStatic();
   }

} 
#endif
