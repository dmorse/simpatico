/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SmpConfigIo.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <mcMd/species/SpeciesMutator.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#ifdef MCMD_LINK
#include <mcMd/links/LinkMaster.h>
#endif
#include <simp/species/Species.h>
#include <util/param/Label.h>
#include <util/param/OptionalLabel.h>
#include <util/misc/FlagSet.h>
#include <util/format/Int.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor. 
   */
   SmpConfigIo::SmpConfigIo(System& system)
    : ConfigIo(system)
   {}

   /* 
   * Destructor.   
   */
   SmpConfigIo::~SmpConfigIo() 
   {}

   /* 
   * Read the configuration file.
   */
   void SmpConfigIo::read(std::istream &in)
   {
      // Make sure Label static buffer is clean on entry
      UTIL_CHECK(Label::isClear());

      using std::endl;

      // Allocate arrays
      int nSpecies = simulation().nSpecies();
      UTIL_CHECK(nSpecies > 0);
      DArray<int> nMoleculeSpecies;
      DArray<int> firstAtomIds;
      nMoleculeSpecies.allocate(nSpecies);
      firstAtomIds.allocate(nSpecies);
      firstAtomIds[0] = 0;

      // Read SPECIES block
      Species* speciesPtr;
      int nAtomSpecies, nAtomMolecule;
      int nAtomTot = 0;
      OptionalLabel speciesBlockLabel("SPECIES"); 
      if (speciesBlockLabel.match(in)) {

         /*
         * If SPECIES block is present, check consistency with data
         * in Species objects, which was read from a parameter file
         * or loaded from an archive upon a restart. Initialize the
         * arrays nMoleculeSpecies and firstAtomIds.
         */

         int iSpecies, nSpeciesIn, iSpeciesIn;
         in >> Label("nSpecies") >> nSpeciesIn;
         UTIL_CHECK(nSpeciesIn > 0);
         UTIL_CHECK(nSpeciesIn == nSpecies);

         Label speciesLabel("species");
         Label nMoleculeLabel("nMolecule");
         for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            in >> speciesLabel>> iSpeciesIn;
            UTIL_CHECK(iSpeciesIn == iSpecies);
            in >> nMoleculeLabel >> nMoleculeSpecies[iSpecies];
            UTIL_CHECK(nMoleculeSpecies[iSpecies] >= 0);
            speciesPtr = &simulation().species(iSpecies);
            UTIL_CHECK(nMoleculeSpecies[iSpecies] <= 
                       speciesPtr->capacity());
            bool match = speciesPtr->matchStructure(in);
            if (!match) {
               UTIL_THROW("Structure mismatch");
            }
            firstAtomIds[iSpecies] = nAtomTot;
            nAtomMolecule = speciesPtr->nAtom();
            nAtomSpecies = nMoleculeSpecies[iSpecies]*nAtomMolecule;
            nAtomTot += nAtomSpecies;
         }

      } else {

         /*
         * If the SPECIES block is absent, use the structures already 
         * defined in the Species objects, and set the nMoleculeSpecies
         * for each species equal to Species::capacity(). Also initialize
         * firstAtomIds array.
         */

         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            speciesPtr = &simulation().species(iSpecies);
            nMoleculeSpecies[iSpecies] = speciesPtr->capacity();
            firstAtomIds[iSpecies+1] = nAtomTot;
            nAtomMolecule = speciesPtr->nAtom();
            nAtomSpecies = nMoleculeSpecies[iSpecies]*nAtomMolecule;
            nAtomTot += nAtomSpecies;
         }
      }

      // Read BOUNDARY block
      in >> Label("BOUNDARY");
      in >> boundary();

      // Read ATOMs block header
      in >> Label("ATOMS");
      OptionalLabel orderedLabel("ordered"); // optional label
      bool isOrdered = orderedLabel.match(in);
      std::string formatString;
      in >> Label("format") >> formatString;
      UTIL_CHECK(formatString.size() > 0);
      int nAtom;
      in >> Label("nAtom") >> nAtom;
      UTIL_CHECK(nAtom == nAtomTot);

      // Parse atom format string
      FlagSet atomFormat("itmpvs");
      atomFormat.setActualOrdered(formatString);
      bool hasAtomIndex = atomFormat.isActive('i');
      bool hasAtomTypeId = atomFormat.isActive('t');
      bool hasAtomContext = atomFormat.isActive('m');
      bool hasAtomPosition = atomFormat.isActive('p');
      bool hasAtomVelocity = atomFormat.isActive('v');
      bool hasAtomShift = atomFormat.isActive('s');
      UTIL_CHECK(hasAtomPosition);
      UTIL_CHECK(isOrdered);

      // TODO: Add ability to read unordered atoms with context.
      // if (!isOrdered) {
      //   UTIL_CHECK(hasAtomContext);
      // }

      // Read all atoms
      Molecule* molPtr;
      Molecule::AtomIterator atomIter;
      int iSpecies, iMol;
      int atomIndex, iAtom, count, atomTypeId;
      int cSpeciesId, cMoleculeId, cAtomId;
      IntVector shift;
      if (isOrdered) {
         count = 0;
         for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            speciesPtr = &simulation().species(iSpecies);
            nAtomMolecule = speciesPtr->nAtom();
            for (iMol = 0; iMol < nMoleculeSpecies[iSpecies]; ++iMol) {

               // Add molecule
               molPtr = &(simulation().getMolecule(iSpecies));
               system().addMolecule(*molPtr);
               if (molPtr != &system().molecule(iSpecies, iMol)) {
                  UTIL_THROW("Molecule index error");
               }

               // Loop over atoms
               iAtom = 0;
               molPtr->begin(atomIter); 
               for ( ; atomIter.notEnd(); ++atomIter) {

                  if (hasAtomIndex) {
                     in >> atomIndex;
                     UTIL_CHECK(atomIndex == count);
                  }
                  if (hasAtomTypeId) {
                     in >> atomTypeId;
                     UTIL_CHECK(atomTypeId 
                                == speciesPtr->atomTypeId(iAtom));
                  } else {
                     atomTypeId = speciesPtr->atomTypeId(iAtom);
                  }
                  atomIter->setTypeId(atomTypeId);
                  if (hasAtomContext) {
                     in >> cSpeciesId;
                     UTIL_CHECK(cSpeciesId == iSpecies);
                     in >> cMoleculeId;
                     UTIL_CHECK(cMoleculeId == iMol);
                     in >> cAtomId;
                     UTIL_CHECK(cAtomId == iAtom);
                  }

                  in >> atomIter->position();
                  if (hasAtomVelocity) {
                     in >> atomIter->velocity();
                  }
                  if (hasAtomShift) {
                     in >> shift;
                     #ifdef MCMD_SHIFT
                     atomIter->shift() = shift;
                     #endif
                  }

                  // Shift atom positions to primary cell
                  #ifdef MCMD_SHIFT
                  boundary().shift(atomIter->position(), atomIter.shift());
                  #else
                  boundary().shift(atomIter->position());
                  #endif

                  iAtom++;
                  count++;
               }  // atom loop
            }  // molecule loop
         }  // species loop
      }  // if (isOrdered)

      #ifdef MCMD_LINK
      if (OptionalLabel("LINKS").match(in)) { 

         // Read Links
         Atom*  atom0Ptr;
         Atom*  atom1Ptr;
         int    iLink, nLink, iAtom, linkType;

         in >> Label("nLink") >> nLink;
         UTIL_CHECK(nLink > 0);
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

      // Make sure static Label buffer is clean on exit
      UTIL_CHECK(Label::isClear());
   }

   /* 
   * Write the configuration file.
   */
   void SmpConfigIo::write(std::ostream &out)
   {
      using std::endl;

      // Write SPECIES block (structure of molecular species)
      out << "SPECIES";
      int nSpecies = simulation().nSpecies();
      out << endl << "nSpecies  " << nSpecies;
      int iSpecies, iMolecule, nMolecule;
      int nAtom, nAtomTot;
      nAtomTot = 0;
      Species* speciesPtr;
      out << endl;
      for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         out << endl << "species " << iSpecies;
         nMolecule = system().nMolecule(iSpecies);
         out << endl << "  nMolecule  " << nMolecule << endl;
         speciesPtr = &simulation().species(iSpecies);
         speciesPtr->writeStructure(out, "  ");
         nAtom = speciesPtr->nAtom();
         nAtomTot += nMolecule*nAtom;
         out << endl;
      }

      // Write BOUNDARY block (boundary dimensions)
      out << endl << "BOUNDARY";
      out << endl << boundary() << endl;

      // Write ATOMS block
      out << endl << "ATOMS";
      out << endl << "ordered";
      out << endl << "format itmpv";
      out << endl << "nAtom " << nAtomTot;
      System::ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator atomIter;
      int atomId, iAtom;
      atomId = 0;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         system().begin(iSpecies, molIter); 
         iMolecule = 0;
         for ( ; molIter.notEnd(); ++molIter) {
            iAtom = 0;
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               out << endl;
               out << atomId << "  ";
               out << atomIter->typeId() << "  ";
               out << iSpecies << " " << iMolecule << " " << iAtom;
               out << endl << atomIter->position();
               out << endl << atomIter->velocity();
               ++iAtom;
               ++atomId;
            }
            ++iMolecule;
         }
      }


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

      UTIL_CHECK(Label::isClear());
   }

} 
