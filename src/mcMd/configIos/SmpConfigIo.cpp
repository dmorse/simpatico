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

      int nSpecies = simulation().nSpecies();
      DArray<int> nMolecules;
      nMolecules.allocate(nSpecies);

      // Read SPECIES block
      in >> Label("SPECIES");
      int nSpeciesIn, iSpeciesIn;
      in >> Label("nSpecies") >> nSpeciesIn;
      UTIL_CHECK(nSpeciesIn > 0);
      UTIL_CHECK(nSpeciesIn == nSpecies);
      Species* speciesPtr;
      Label speciesLabel("species");
      Label nMoleculeLabel("nMolecule");
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         in >> Label("species") >> iSpeciesIn;
         UTIL_CHECK(iSpeciesIn == iSpecies);
         in >> Label("nMolecule") >> nMolecules[iSpecies];
         speciesPtr = &simulation().species(iSpecies);
         UTIL_CHECK( nMolecules[iSpecies] <= speciesPtr->capacity() );
         //speciesPtr->checkStructure(in);
      }
      
      // Read BOUNDARY block
      in >> Label("BOUNDARY");
      in >> boundary();


      // Read ATOM block
      in >> Label("ATOM");
      Molecule::AtomIterator atomIter;
      Molecule* molPtr;
      Label moleculeLabel("molecule");
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         in >> Label("species") >> iSpeciesIn;
         if (iSpeciesIn != iSpecies) {
            UTIL_THROW("Error: iSpeciesIn != iSpecies");
         }
         for (int iMol = 0; iMol < nMolecules[iSpecies]; ++iMol) {
            molPtr = &(simulation().getMolecule(iSpecies));
            system().addMolecule(*molPtr);
            if (molPtr != &system().molecule(iSpecies, iMol)) {
               UTIL_THROW("Molecule index error");
            }
            for (molPtr->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               // readAtom
            }

         }
      }

   }

   /* 
   * Write the configuration file.
   */
   void SmpConfigIo::write(std::ostream &out)
   {
      using std::endl;

      // Write Boundary dimensions
      out << "BOUNDARY" << endl << endl;
      out << boundary() << endl;

      // Species atom info
      out << endl << "SPECIES" << endl;
      int nSpecies = simulation().nSpecies();
      out << endl << "nSpecies  " << nSpecies;
      int iMolecule, nMolecule;
      int nAtom, nAtomTot;
      nAtomTot = 0;
      Species* speciesPtr;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         out << endl << "species " << iSpecies;
         nMolecule = system().nMolecule(iSpecies);
         out << endl << "  nMolecule  " << nMolecule;
         speciesPtr = &simulation().species(iSpecies);
         speciesPtr->writeStructure(out, "  ");
         nAtom = speciesPtr->nAtom();
         nAtomTot += nMolecule*nAtom;
      }

      // Write ATOM information
      out << endl << endl << "ATOM" << endl;
      out << endl << "ordered";
      out << endl << "atomFormat imtpv";
      System::ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator   atomIter;
      int atomId, iAtom;
      atomId = 0;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         out << endl;
         speciesPtr = &simulation().species(iSpecies);
         system().begin(iSpecies, molIter); 
         iMolecule = 0;
         for ( ; molIter.notEnd(); ++molIter) {
            iAtom = 0;
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               out << endl;
               out << atomId << " ";
               out << iSpecies << " " << iMolecule
                   << " " << iAtom << " ";
               out << atomIter->typeId() << " ";
               out << atomIter->position() << " ";
               out << atomIter->velocity() << " ";
               ++iAtom;
               ++atomId;
            }
            ++iMolecule;
         }
      }
   }

} 
