/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpcConfigIo.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <simp/species/Species.h>
#include <mcMd/species/SpeciesMutator.h>
#include <mcMd/chemistry/Group.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#ifdef SIMP_TETHER
#include <mcMd/tethers/TetherMaster.h>
#endif
#ifdef MCMD_LINK
#include <mcMd/links/LinkMaster.h>
#endif
#include <util/param/Label.h>
#include <util/misc/FlagSet.h>
#include <util/format/Int.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor. 
   */
   SpcConfigIo::SpcConfigIo(System& system)
    : ConfigIo(system)
   {}

   /* 
   * Destructor.   
   */
   SpcConfigIo::~SpcConfigIo() 
   {}

   /* 
   * Read the configuration file.
   */
   void SpcConfigIo::read(std::istream &in)
   {

      // System boundary.
      in >> Label("BOUNDARY");
      in >> boundary();

      // Molecular species
      in >> Label("SPECIES");
      #ifdef SIMP_BOND
      bool hasBonds = (bool)simulation().nBondType();
      if (hasBonds) {
         Label label("hasBonds", false);
         in >> label;
         hasBonds = label.isClear();
      }
      #endif
      #ifdef SIMP_ANGLE
      bool hasAngles = (bool)simulation().nAngleType();
      if (hasAngles) {
         Label label("hasAngles", false);
         in >> label;
         hasAngles = label.isClear();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      bool hasDihedrals = (bool)simulation().nDihedralType();
      if (hasDihedrals) {
         Label label("hasDihedrals", false);
         in >> label;
         hasDihedrals = label.isClear();
      }
      #endif
      int nSpecies = simulation().nSpecies();
      int nSpeciesIn;
      in >> Label("nSpecies") >> nSpeciesIn;
      UTIL_CHECK(nSpeciesIn == nSpecies);
      Species* speciesPtr;
      Molecule* molPtr;
      int iSpeciesIn, nMolecule, nAtomIn;
      int nAtomTot = 0;
      for (int iSpecies = 0 ; iSpecies < nSpecies; ++iSpecies) {
         in >> iSpeciesIn;
         UTIL_CHECK(iSpeciesIn == iSpecies);
         speciesPtr = &simulation().species(iSpecies);
         in >> nMolecule;
         UTIL_CHECK(nMolecule <= speciesPtr->capacity()) 
         in >> nAtomIn;
         UTIL_CHECK(nAtomIn == speciesPtr->nAtom())
         nAtomTot += nMolecule*nAtomIn;
         #ifdef SIMP_BOND
         if (hasBonds) {
            int nBondIn;
            in >> nBondIn;
            UTIL_CHECK(nBondIn == speciesPtr->nBond());
         }
         #endif
         #ifdef SIMP_ANGLE
         if (hasAngles) {
            int nAngleIn;
            in >> nAngleIn;
            UTIL_CHECK(nAngleIn == speciesPtr->nAngle());
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (hasDihedrals) {
            int nDihedralIn;
            in >> nDihedralIn;
            UTIL_CHECK(nDihedralIn == speciesPtr->nDihedral());
         }
         #endif
         // Add all molecules of this species
         for (int iMol = 0; iMol < nMolecule; ++iMol) {
            molPtr = &(simulation().getMolecule(iSpecies));
            system().addMolecule(*molPtr);
            if (molPtr != &system().molecule(iSpecies, iMol)) {
               UTIL_THROW("Molecule index error");
            }
         }
      }

      // Atoms block
      in >> Label("ATOMS");
      Label orderedLabel("ordered", false);
      in >> orderedLabel;
      bool isOrdered = orderedLabel.isClear();
      std::string format;
      in >> Label("format") >> format;
      FlagSet flags("imtpvsc");
      flags.setActualOrdered(format);
      UTIL_CHECK(flags.isActive('i'));
      UTIL_CHECK(flags.isActive('m'));
      UTIL_CHECK(flags.isActive('t'));
      UTIL_CHECK(flags.isActive('p'));
      bool hasVelocity = flags.isActive('v');
      //bool hasShifts = flags.isActive('s');
      in >> Label("nAtom") >> nAtomIn;
      UTIL_CHECK(nAtomIn == nAtomTot);
      
      if (isOrdered) {   
         System::MoleculeIterator molIter;
         Molecule::AtomIterator atomIter;
         int atomId, atomIdIn, iMol, iMolIn, iAtom, iAtomIn, typeId;
         atomId = 0;
         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            speciesPtr = &simulation().species(iSpecies);
            iMol = 0;
            system().begin(iSpecies, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               iAtom = 0;
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                  in >> atomIdIn;
                  UTIL_CHECK(atomIdIn == atomId);
                  in >> iSpeciesIn;
                  UTIL_CHECK(iSpeciesIn == iSpecies);
                  in >> iMolIn;
                  UTIL_CHECK(iMolIn == iMol);
                  in >> iAtomIn;
                  UTIL_CHECK(iAtomIn == iAtom);
                  in >> typeId;
                  atomIter->setTypeId(typeId);
                  in >> atomIter->position();
                  if (hasVelocity) {
                     in >> atomIter->velocity();
                  }
                  ++atomId;
                  ++iAtom;
               }
               ++iMol; 
            }
         }
      }

      #ifdef SIMP_TETHER
      {  // Scope block for variables 

         // Read Tethers
         Vector anchor;
         Atom*  atomPtr;
         int    iTether, nTether, iAtom;
         in >> Label("TETHERS");
         in >> Label("nTether") >> nTether;
         for (iTether = 0; iTether < nTether; ++iTether) {
            in >> iSpecies >> iMol >> iAtom >> anchor;
            if (iSpecies > = simulation().nSpecies()) {
               UTIL_THROW("Invalid species index for tether");
            }
            if (iMol > = system().nMolecule(iSpecies)) {
               UTIL_THROW("Invalid molecule index for tether");
            }
            molPtr = &system().molecule(iSpecies, iMol);
            if (iAtom > = molPtr->nAtom()) {
               Log::file() << iAtom << "   " << molPtr->nAtom() << std::endl;
               UTIL_THROW("Invalid atom index for tether");
            }
            atomPtr = &molPtr->atom(iAtom);
            system().tetherMaster().addTether(*atomPtr, anchor);
         }

      }
      #endif

      #ifdef MCMD_LINK
      { // Scope block for variables

         // Read Links
         Atom*  atom0Ptr;
         Atom*  atom1Ptr;
         int    iLink, nLink, iAtom, linkType;

         in >> Label("LINKS");
         in >> Label("nLink") >> nLink;
         for (iLink = 0; iLink < nLink; ++iLink) {

            // Read atom 0
            in >> iSpecies >> iMol >> iAtom;
            if (iSpecies < 0 || iSpecies > = simulation().nSpecies()) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iSpecies0 index in link");
            }
            if (iMol < 0 || iMol > = system().nMolecule(iSpecies)) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iMol0 index in link");
            }
            molPtr = &system().molecule(iSpecies, iMol);
            if (iAtom < 0 || iAtom > = molPtr->nAtom()) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iAtom0 index for link");
            }
            atom0Ptr = &molPtr->atom(iAtom);

            // Read atom 1
            in >> iSpecies >> iMol >> iAtom;
            if (iSpecies < 0 || iSpecies > = simulation().nSpecies()) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iSpecies1 index in link");
            }
            if (iMol < 0 || iMol > = system().nMolecule(iSpecies)) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iMol1 index in link");
            }
            molPtr = &system().molecule(iSpecies, iMol);
            if (iAtom < 0 || iAtom > = molPtr->nAtom()) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iAtom1 index for link");
            }
            atom1Ptr = &molPtr->atom(iAtom);

            in >> linkType;
            system().linkMaster().addLink(*atom0Ptr, *atom1Ptr, linkType);

         }

      }
      #endif

   }

   /* 
   * Write the configuration file.
   */
   void SpcConfigIo::write(std::ostream &out)
   {
      using std::endl;

      // Write Boundary dimensions
      out << "BOUNDARY\n";
      out << boundary() << endl;

      // Write species information
      out << endl << "SPECIES\n";
      #ifdef SIMP_BOND
      bool hasBonds = (bool)simulation().nBondType();
      int  nBondTot = 0;
      if (hasBonds) {
         out << "hasBonds\n";
      }
      #endif
      #ifdef SIMP_ANGLE
      bool hasAngles = (bool)simulation().nAngleType();
      int  nAngleTot = 0;
      if (hasAngles) {
         out << "hasAngles\n";
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      bool hasDihedrals = (bool)simulation().nDihedralType();
      int  nDihedralTot = 0;
      if (hasDihedrals) {
         out << "hasDihedrals\n";
      }
      #endif
      Species* speciesPtr;
      int iSpecies, nSpecies, nMolecule, nAtomSpecies, nAtomTot;
      nSpecies = simulation().nSpecies();
      out << "nSpecies  " << nSpecies << "\n";
      nAtomTot = 0;
      for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         out << iSpecies;
         out << "  " << speciesPtr->capacity();
         nAtomSpecies = speciesPtr->nAtom();
         out << "  " << nAtomSpecies;
         nMolecule = system().nMolecule(iSpecies);
         nAtomTot += nAtomSpecies*nMolecule;
         #ifdef SIMP_BOND
         if (hasBonds) {
            out << "  " << speciesPtr->nBond();
            nBondTot += nMolecule*(speciesPtr->nBond());
         }
         #endif
         #ifdef SIMP_ANGLE
         if (hasAngles) {
            out << "  " << speciesPtr->nAngle();
            nAngleTot += nMolecule*(speciesPtr->nAngle());
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (hasDihedrals) {
            out << "  " << speciesPtr->nDihedral();
            nDihedralTot += nMolecule*(speciesPtr->nDihedral());
         }
         #endif
         out << "\n";
      }
      out << "\n";

      // Write atomic positions
      out << "ATOMS\n";
      out << "ordered\n";
      out << "format ";
      std::string format = "imtp";
      out << format << "\n";
      out << "nAtom " << nAtomTot << "\n";

      System::ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator atomIter;
      int iMol, iAtom;
      int atomId = 0;
      for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         iMol = 0;
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            iAtom = 0;
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               out << atomId << "  ";
               out << iSpecies << "  ";
               out << iMol << "  ";
               out << iAtom << "  ";
               out << atomIter->typeId() << "  ";
               out << atomIter->position() << "  ";
               out << "\n";
               ++iAtom;
               ++atomId;
            }
            ++iMol;
         }
      }
      out << "\n";

      #ifdef SIMP_BOND
      if (hasBonds) {
         out << "BONDS\n";
         out << "nBond  " << nBondTot << "\n";
         writeGroups<2>(out);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAngles) {
         out << "ANGLES\n";
         out << "nAngle  " << nAngleTot << "\n";
         writeGroups<3>(out);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedrals) {
         out << "nDihedral  " << nDihedralTot << "\n";
         out << "DIHEDRALS\n";
         writeGroups<4>(out);
      }
      #endif

      #if 0
      #ifdef SIMP_TETHER
      { // Scope for local variables

         // Write Tethers
         Tether*   tetherPtr;
         Atom*     atomPtr;
         Molecule* molPtr;
         int       iTether, nTether, iAtom;
         out << std::endl;
         out << "TETHERS" << endl << endl;
         nTether = system().tetherMaster().nTether();
         out << Label("nTether") << nTether << std::endl;
         for (iTether = 0; iTether < nTether; ++iTether) {
            tetherPtr = &(system().tetherMaster().tether(iTether));
            atomPtr = &tetherPtr->atom();
            molPtr = &atomPtr->molecule();
            iAtom = atomPtr->indexInMolecule();
            iMol = system().moleculeId(*molPtr);
            iSpecies = molPtr->species().id();
            out << Int(iSpecies,5) << Int(iMol,9) << Int(iAtom,6) 
                << tetherPtr->anchor() << std::endl;
         }

      }
      #endif
      #endif

      #if 0
      #ifdef MCMD_LINK
      { // Scope for local variables

         // Write Links
         Link*      linkPtr;
         Atom*      atomPtr;
         Molecule*  molPtr;
         int        iLink, nLink, iAtom;
         out << std::endl;
         out << "LINKS" << endl << endl;
         nLink = system().linkMaster().nLink();
         out << Label("nLink") << nLink << std::endl;
         for (iLink = 0; iLink < nLink; ++iLink) {
            linkPtr = &(system().linkMaster().link(iLink));

            // Output species, molecule, atom ids for atom 0
            atomPtr = &(linkPtr->atom0());
            molPtr = &atomPtr->molecule();
            iAtom = atomPtr->indexInMolecule();
            iMol = system().moleculeId(*molPtr);
            iSpecies = molPtr->species().id();
            out << Int(iSpecies,8) << Int(iMol,8) << Int(iAtom,8);
            out << "   ";

            // Output species, molecule, atom ids for atom 1
            atomPtr = &(linkPtr->atom1());
            molPtr = &atomPtr->molecule();
            iAtom = atomPtr->indexInMolecule();
            iMol = system().moleculeId(*molPtr);
            iSpecies = molPtr->species().id();
            out << Int(iSpecies,8) << Int(iMol,8) << Int(iAtom,8);
            out << "   ";

            out << Int(linkPtr->typeId(),8) << std::endl;
         }

      }
      #endif
      #endif

   }

   /* 
   * Write the configuration file.
   */
   template <int N>
   void SpcConfigIo::writeGroups(std::ostream &out)
   {
      System::ConstMoleculeIterator molIter;
      ConstArrayIterator< Group<N> > groupIter;
      int i;
      int nSpecies = simulation().nSpecies();
      int groupId = 0;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(groupIter); groupIter.notEnd(); ++groupIter) {
               out << "  " << groupId;
               out << "  " << groupIter->typeId();
               for (i = 0; i < N; ++i) {
                  out << "  " << groupIter->atom(i).id();
               }
               out << "\n";
               ++groupId;
            }
         }
      }
      out << "\n";
   }

} 
