/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpcConfigIo.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <mcMd/species/Species.h>
#include <mcMd/species/SpeciesMutator.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#ifdef INTER_TETHER
#include <mcMd/tethers/TetherMaster.h>
#endif
#ifdef MCMD_LINK
#include <mcMd/links/LinkMaster.h>
#endif
#include <util/param/Label.h>
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
      in >> boundary() ;

      // Molecular species
      Species* speciesPtr;
      Molecule* molPtr;
      int iSpecies, iSpeciesIn, nMolecule, nAtom, nAtomTot;
      int nSpecies = simulation().nSpecies();
      bool isRequired = false;
      bool isMolecular = false;
      in >> Label("SPECIES", isRequired);
      if (Label::isClear()) {
         isMolecular = true;
         nAtomTot = 0;
         iSpecies = 0;
         for ( ; iSpecies < nSpecies; ++iSpecies) {
            in >> iSpeciesIn;
            if (iSpeciesIn != iSpecies) {
               UTIL_THROW("Error: iSpeciesIn != iSpecies");
            }
            speciesPtr = &simulation().species(iSpecies);
            in >> nMolecule;
            in >> nAtom;
            if (nMolecule > speciesPtr->capacity()) {
               UTIL_THROW("Error: nMolecule > species.capacity()");
            }
            if (nAtom != speciesPtr->nAtom()) {
               UTIL_THROW("Error: nAtom != species.nAtom() ");
            }
            // Add all molecules of this species
            for (iMol = 0; iMol < nMolecule; ++iMol) {
               molPtr = &(simulation().getMolecule(iSpecies));
               system().addMolecule(*molPtr);
               if (molPtr != &system().molecule(iSpecies, iMol)) {
                  UTIL_THROW("Molecule index error");
               }
            }
            nAtomTot += nMolecule*nAtom;
         }
      }

      in >> Label("ATOMS");
      int nAtom, atomId, iMol; iAtom;
      in >> Label("nAtom") >> nAtom;
      Molecule::AtomIterator atomIter;
      for (iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         in >> atomId;
         in >> iSpeciesIn;
         if (iSpeciesIn != iSpecies) {
            UTIL_THROW("Error: iSpeciesIn != iSpecies");
         }
         in >> Label("nMolecule") >> nMolecule;
         for (iMol = 0; iMol < nMolecule; ++iMol) {
            in >> moleculeLabel >> iMolIn;
            molPtr = &(simulation().getMolecule(iSpecies));
            system().addMolecule(*molPtr);

            if (iMolIn != iMol) {
               UTIL_THROW("Error: iMolIn != iMol");
            }

            // If mutable, read molecule state id
            if (speciesPtr->isMutable()) {
               speciesPtr->mutator().readMoleculeState(in, *molPtr);
            }

            // Read positions, shift into primary cell
            for (molPtr->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               readAtom(in, *atomIter);
            }

         }
      }

      #ifdef INTER_TETHER
      {  // Scope block for variables 

         // Read Tethers
         Vector anchor;
         Atom*  atomPtr;
         int    iTether, nTether, iAtom;
         in >> Label("TETHERS");
         in >> Label("nTether") >> nTether;
         for (iTether = 0; iTether < nTether; ++iTether) {
            in >> iSpecies >> iMol >> iAtom >> anchor;
            if (iSpecies >= simulation().nSpecies()) {
               UTIL_THROW("Invalid species index for tether");
            }
            if (iMol >= system().nMolecule(iSpecies)) {
               UTIL_THROW("Invalid molecule index for tether");
            }
            molPtr  = &system().molecule(iSpecies, iMol);
            if (iAtom >= molPtr->nAtom()) {
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
            if (iSpecies < 0 || iSpecies >= simulation().nSpecies()) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iSpecies0 index in link");
            }
            if (iMol < 0 || iMol >= system().nMolecule(iSpecies)) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iMol0 index in link");
            }
            molPtr  = &system().molecule(iSpecies, iMol);
            if (iAtom < 0 || iAtom >= molPtr->nAtom()) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iAtom0 index for link");
            }
            atom0Ptr = &molPtr->atom(iAtom);

            // Read atom 1
            in >> iSpecies >> iMol >> iAtom;
            if (iSpecies < 0 || iSpecies >= simulation().nSpecies()) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iSpecies1 index in link");
            }
            if (iMol < 0 || iMol >= system().nMolecule(iSpecies)) {
               Log::file() << "iLink = " << iLink << std::endl;
               UTIL_THROW("Invalid iMol1 index in link");
            }
            molPtr  = &system().molecule(iSpecies, iMol);
            if (iAtom < 0 || iAtom >= molPtr->nAtom()) {
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
      out << "BOUNDARY" << endl << endl;
      out << boundary() << endl;

      // Write species information
      Species* speciesPtr;
      int iSpecies, nSpecies;
      out << endl << "SPECIES" << endl;
      nSpecies = simulation().nSpecies;
      out << "nSpecies   " << iSpecies << endl;
      for (iSpecies = 0; iSpecies < simulation().nSpecies(); ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         out << iSpecies << "  ";
         out << speciesPtr->capacity() << "  ";
         out << speciesPtr->nAtom() << endl;
      }

      // Write atomic positions
      System::ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator   atomIter;
      Species* speciesPtr;
      int iSpecies, iMolecule, iAtom, atomId;
      out << endl << "ATOMS" << endl << endl;
      atomId = 0;
      iSpecies = 0;
      for ( ; iSpecies < simulation().nSpecies(); ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         iMolecule = 0;
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            // if (speciesPtr->isMutable()) {
            //   speciesPtr->mutator().writeMoleculeState(out, *molIter);
            // }
            iAtom = 0;
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               writeAtom(out, *atomIter);
               out << atomId << "  ";
               out << iSpecies << "  ";
               out << iMolecule << "  ";
               out << iAtom << "  ";
               out << atomIter->typeId() << "  ";
               out << atomIter->position() << "  ";
               out << endl;
               ++iAtom;
            }
            ++iMolecule;
         }
      }

      #if 0
      #ifdef INTER_TETHER
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
            atomPtr   = &tetherPtr->atom();
            molPtr    = &atomPtr->molecule();
            iAtom     = atomPtr->indexInMolecule();
            iMolecule = system().moleculeId(*molPtr);
            iSpecies  = molPtr->species().id();
            out << Int(iSpecies,5) << Int(iMolecule,9) << Int(iAtom,6) 
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
            linkPtr  = &(system().linkMaster().link(iLink));

            // Output species, molecule, atom ids for atom 0
            atomPtr   = &(linkPtr->atom0());
            molPtr    = &atomPtr->molecule();
            iAtom     = atomPtr->indexInMolecule();
            iMolecule = system().moleculeId(*molPtr);
            iSpecies  = molPtr->species().id();
            out << Int(iSpecies,8) << Int(iMolecule,8) << Int(iAtom,8);
            out << "   ";

            // Output species, molecule, atom ids for atom 1
            atomPtr   = &(linkPtr->atom1());
            molPtr    = &atomPtr->molecule();
            iAtom     = atomPtr->indexInMolecule();
            iMolecule = system().moleculeId(*molPtr);
            iSpecies  = molPtr->species().id();
            out << Int(iSpecies,8) << Int(iMolecule,8) << Int(iAtom,8);
            out << "   ";

            out << Int(linkPtr->typeId(),8) << std::endl;
         }

      }
      #endif
      #endif

   }

} 
