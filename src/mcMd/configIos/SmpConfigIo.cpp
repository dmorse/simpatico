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
#ifdef SIMP_TETHER
#include <mcMd/tethers/TetherMaster.h>
#endif
#ifdef MCMD_LINK
#include <mcMd/links/LinkMaster.h>
#endif
#include <simp/species/Species.h>
#include <util/param/Label.h>
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
      DArray<int> firstAtomIndex;
      nMoleculeSpecies.allocate(nSpecies);
      firstAtomIndex.allocate(nSpecies);
      firstAtomIndex[0] = 0;

      // Read SPECIES block
      Species* speciesPtr;
      int nAtomSpecies, nAtomMolecule;
      int nAtomTot = 0;
      in >> Label("SPECIES");
      if (Label::isClear()) {
         int nSpeciesIn, iSpeciesIn;
         in >> Label("nSpecies") >> nSpeciesIn;
         UTIL_CHECK(nSpeciesIn > 0);
         UTIL_CHECK(nSpeciesIn == nSpecies);
         Label speciesLabel("species");
         Label nMoleculeLabel("nMolecule");
         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            in >> Label("species") >> iSpeciesIn;
            UTIL_CHECK(iSpeciesIn == iSpecies);
            in >> Label("nMolecule") >> nMoleculeSpecies[iSpecies];
            UTIL_CHECK(nMoleculeSpecies[iSpecies] >= 0);
            speciesPtr = &simulation().species(iSpecies);
            UTIL_CHECK(nMoleculeSpecies[iSpecies] <= 
                       speciesPtr->capacity());
            bool match = speciesPtr->matchStructure(in);
            if (!match) {
               UTIL_THROW("Structure mismatch");
            }
            nAtomMolecule = speciesPtr->nAtom();
            nAtomSpecies = nMoleculeSpecies[iSpecies]*nAtomMolecule;
            nAtomTot += nAtomSpecies;
            if (iSpecies < nSpecies - 1) {
               firstAtomIndex[iSpecies+1] = firstAtomIndex[iSpecies] 
                                          + nAtomSpecies;
            }
         }
      } else {
         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            speciesPtr = &simulation().species(iSpecies);
            nMoleculeSpecies[iSpecies] = speciesPtr->capacity();
            nAtomMolecule = speciesPtr->nAtom();
            nAtomSpecies = nMoleculeSpecies[iSpecies]*nAtomMolecule;
            nAtomTot += nAtomSpecies;
            if (iSpecies < nSpecies - 1) {
               firstAtomIndex[iSpecies+1] = firstAtomIndex[iSpecies] 
                                          + nAtomSpecies;
            }
         }
      }

      // Read BOUNDARY block
      in >> Label("BOUNDARY");
      in >> boundary();

      // Read ATOM block header
      in >> Label("ATOM");
      in >> Label("ordered");
      bool isOrdered = Label::isClear();
      std::string formatString;
      in >> Label("format");
      in >> formatString;
      UTIL_CHECK(formatString.size() > 0);
      int nAtom;
      in >> Label("nAtom") >> nAtom;
      UTIL_CHECK(nAtom == nAtomTot);

      // Parse atom format string
      FlagSet atomFormat("imtpvs");
      atomFormat.setActualOrdered(formatString);
      bool hasAtomIndex = atomFormat.isActive('i');
      bool hasAtomContext = atomFormat.isActive('m');
      bool hasAtomTypeId = atomFormat.isActive('t');
      bool hasAtomPosition = atomFormat.isActive('p');
      bool hasAtomVelocity = atomFormat.isActive('v');
      bool hasAtomShift = atomFormat.isActive('s');
      UTIL_CHECK(hasAtomPosition);
      if (!isOrdered) {
         UTIL_CHECK(hasAtomContext);
      }

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

                  // Log::file() << endl;
                  if (hasAtomIndex) {
                     in >> atomIndex;
                     UTIL_CHECK(atomIndex == count);
                     // Log::file() << atomIndex << " ";
                  }
                  if (hasAtomContext) {
                     in >> cSpeciesId;
                     UTIL_CHECK(cSpeciesId == iSpecies);
                     // Log::file() << cSpeciesId << " ";
                     in >> cMoleculeId;
                     UTIL_CHECK(cMoleculeId == iMol);
                     // Log::file() << cMoleculeId << " ";
                     in >> cAtomId;
                     UTIL_CHECK(cAtomId == iAtom);
                     // Log::file() << cAtomId << " ";
                  }
                  if (hasAtomTypeId) {
                     in >> atomTypeId;
                     UTIL_CHECK(atomTypeId 
                                == speciesPtr->atomTypeId(iAtom));
                     // Log::file() << atomTypeId << " ";
                  } else {
                     atomTypeId = speciesPtr->atomTypeId(iAtom);
                  }
                  atomIter->setTypeId(atomTypeId);

                  in >> atomIter->position();
                  // Log::file() << endl << atomIter->position(); 
                  if (hasAtomVelocity) {
                     in >> atomIter->velocity();
                     // Log::file() << endl << atomIter->velocity();
                  }
                  if (hasAtomShift) {
                     in >> shift;
                     #ifdef MCMD_SHIFT
                     atomIter->shift() = shift;
                     #endif
                  }

                  // Shift atom positions to primary cell
                  #ifdef MCMD_SHIFT
                  in >> atomIter->shift();
                  boundary().shift(atomIter->position(), atomIter.shift());
                  #else
                  boundary().shift(atomIter->position());
                  #endif

                  iAtom++;
                  count++;
               }   // atom loop
            }   // molecule loop
         }   // species loop
      }   // if (isOrdered)

      #if 0
      // Add all molecules
      int nMolecule, iMol;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         nMolecule = nMoleculeSpecies[iSpecies];
         for (int iMol = 0; iMol < nMolecule; ++iMol) {
            molPtr = &(simulation().getMolecule(iSpecies));
            system().addMolecule(*molPtr);
            if (molPtr != &system().molecule(iSpecies, iMol)) {
               UTIL_THROW("Molecule index error");
            }
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

      // Species atom info
      out << "SPECIES";
      int nSpecies = simulation().nSpecies();
      out << endl << "nSpecies  " << nSpecies;
      int iMolecule, nMolecule;
      int nAtom, nAtomTot;
      nAtomTot = 0;
      Species* speciesPtr;
      out << endl;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         out << endl << "species " << iSpecies;
         nMolecule = system().nMolecule(iSpecies);
         out << endl << "  nMolecule  " << nMolecule;
         speciesPtr = &simulation().species(iSpecies);
         speciesPtr->writeStructure(out, "  ");
         nAtom = speciesPtr->nAtom();
         nAtomTot += nMolecule*nAtom;
         out << endl;
      }

      // Write Boundary dimensions
      out << endl << "BOUNDARY";
      out << endl << boundary() << endl;

      // Write ATOM information
      out << endl << "ATOM";
      out << endl << "ordered";
      out << endl << "format imtpv";
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
               out << iSpecies << " " << iMolecule << " " << iAtom << "  ";
               out << atomIter->typeId();
               out << endl << atomIter->position();
               out << endl << atomIter->velocity();
               ++iAtom;
               ++atomId;
            }
            ++iMolecule;
         }
      }

      UTIL_CHECK(Label::isClear());
   }

} 
