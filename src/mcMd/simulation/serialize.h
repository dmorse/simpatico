#ifndef MCMD_SYSTEM_SERIALIZE_H
#define MCMD_SYSTEM_SERIALIZE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Serialize a System configuration to / from an archive.
   */
   template <class Archive>
   void System::serialize(Archive& ar, const unsigned int version)
   {
      ar & boundary();

      Species* speciesPtr;
      int nSpecies = simulation().nSpecies();

      if (Archive::is_saving()) {

         System::MoleculeIterator molIter;
         Molecule::AtomIterator   atomIter;
         int nMoleculeOut;

         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            ar & iSpecies;
            nMoleculeOut = nMolecule(iSpecies);
            ar & nMoleculeOut;
            speciesPtr = &simulation().species(iSpecies);
            begin(iSpecies, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter); 
               for ( ; atomIter.notEnd(); ++atomIter) {
                  #ifdef MCMD_SHIFT
                  boundary().shift(atomIter->position(), atomIter->shift());
                  ar & atomIter->shift();
                  #else
                  boundary().shift(atomIter->position());
                  #endif
                  ar & atomIter->position();
                  ar & atomIter->velocity();
               }
            }
         }

      }

      if (Archive::is_loading()) {

         Molecule* molPtr;
         Molecule::AtomIterator atomIter;
         int iSpeciesIn, nMoleculeIn;

         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            speciesPtr = &simulation().species(iSpecies);
            ar & iSpeciesIn;
            if (iSpeciesIn != iSpecies) {
               UTIL_THROW("Error: iSpeciesIn != iSpecies");
            }
            ar & nMoleculeIn;
            for (int iMol = 0; iMol < nMoleculeIn; ++iMol) {
               molPtr = &(speciesPtr->reservoir().pop());
               addMolecule(*molPtr);
               if (molPtr != &molecule(iSpecies, iMol)) {
                  UTIL_THROW("Molecule index error");
               }
               molPtr->begin(atomIter); 
               for ( ; atomIter.notEnd(); ++atomIter) {
                  ar & atomIter->position();
                  ar & atomIter->velocity();
                  #ifdef MCMD_SHIFT
                  ar & atomIter->shift();
                  boundary().shift(atomIter->position(), atom.shift());
                  #else
                  boundary().shift(atomIter->position());
                  #endif
               }
            }
         }
      }

   }

} 
#endif
