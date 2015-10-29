/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Linear.h"

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/neighbor/CellList.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/chemistry/Molecule.h>

#include <sstream>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   Linear::Linear() 
   : Species(),
     #ifdef INTER_ANGLE
     hasAngles_(0),
     #endif
     #ifdef INTER_DIHEDRAL
     hasDihedrals_(0),
     #endif
     maxPlacementAttempts_(200)
   {}

   /*
   * Destructor.
   */
   Linear::~Linear() 
   {}

   /*
   * Build the chemical structure of a flexible chain molecule.
   */
   void Linear::buildLinear() 
   {

      // Preconditions
      if (nAtom() < 2) {
         UTIL_THROW("nAtom < 2");
      }
      if (nBond() != nAtom() - 1) {
         UTIL_THROW("nBond != nAtom - 1");
      }
      #ifdef INTER_ANGLE
      if (hasAngles_) {
         if (nAngle() != nAtom() - 2) {
            UTIL_THROW("nAngle != nAtom - 2");
         }
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedrals_) {
         if (nAtom() > 3) {
            if (nDihedral() != nAtom() - 3) {
               UTIL_THROW("nDihedral != nAtom - 3");
            }
         } else {
            if (nDihedral() != 0) {
               UTIL_THROW("nDihedral != 0 when nAtom < 4");
            }
         }
      }
      #endif

      int i; 

      allocate();

      // Set Atom Type Ids
      for (i = 0; i < nAtom(); ++i) {
         setAtomType(i, calculateAtomTypeId(i));
      }

      // Build Bonds
      for (i = 0; i < nBond(); ++i) {
         makeBond(i, i, i+1, calculateBondTypeId(i));
      }

      #ifdef INTER_ANGLE
      // Build Angles 
      if (hasAngles_) {
         for (i = 0; i < nAngle(); ++i) {
            makeAngle(i, i, i+1, i+2, calculateAngleTypeId(i));
         }
      }
      #endif

      #ifdef INTER_DIHEDRAL
      // Build Dihedrals.
      if (hasDihedrals_) {
         for (i = 0; i < nDihedral(); ++i) {
            makeDihedral(i, i, i+1, i+2, i+3, calculateDihedralTypeId(i));
         }
      }
      #endif

   }

   /**
   * Recursive function to try to place an atom.
   */
   bool Linear::tryPlaceAtom(Molecule& molecule, int atomId,
      DArray<double> exclusionRadius, System& system, CellList &cellList,
      BondPotential *bondPotentialPtr, const Boundary &boundary)
   {
      Atom& lastAtom = molecule.atom(atomId);
      Atom& thisAtom = molecule.atom(++atomId);
      Random& random = system.simulation().random();

      bool hasBeenPlaced = false;

      for (int iAttempt = 0; iAttempt < maxPlacementAttempts_; iAttempt++)
      {
         // draw a random bond vector
         int beta = 1;
         Vector v;
         random.unitVector(v);
         v *= bondPotentialPtr->randomBondLength(&random, beta,
            calculateBondTypeId(lastAtom.indexInMolecule()));

         Vector newPos;
         newPos = lastAtom.position();
         newPos += v;
         // shift into simulation cell
         boundary.shift(newPos);

         // check if the atom can be placed at the new position
         CellList::NeighborArray neighbors;
         cellList.getNeighbors(newPos, neighbors);
         int nNeighbor = neighbors.size();
         bool canBePlaced = true;
         for (int j = 0; j < nNeighbor; ++j) {
            Atom *jAtomPtr = neighbors[j];
         
            double r = sqrt(boundary.distanceSq(
               jAtomPtr->position(), newPos));
            if (r < (exclusionRadius[thisAtom.typeId()] +
               exclusionRadius[jAtomPtr->typeId()])) {
               canBePlaced = false;
               break;
            }
         } 
         if (canBePlaced) {
            // place the particle
            thisAtom.position() = newPos;

            // add to cell list
            cellList.addAtom(thisAtom);

            // are we add the end of the chain?
            if (atomId == molecule.nAtom()-1) 
              return true;
            
            // recursion step
            if (! tryPlaceAtom(molecule, atomId, exclusionRadius, system,
               cellList, bondPotentialPtr, boundary) ) {
               // If the next monomer cannot be inserted, delete this monomer
               // again
               cellList.deleteAtom(thisAtom);
            } else {
               hasBeenPlaced = true;
               break;
            }
         } 
      }

      return hasBeenPlaced;
   }

   /**
   * Generate random molecules
   */
   void Linear::generateMolecules(int nMolecule, 
      DArray<double> exclusionRadius, System& system,
      BondPotential *bondPotentialPtr, const Boundary &boundary)
   {
      int iMol;

      // Set up a cell list with twice the maxium exclusion radius as the
      // cell size
      double maxExclusionRadius = 0.0;
      for (int iType = 0; iType < system.simulation().nAtomType(); iType++) {
         if (exclusionRadius[iType] > maxExclusionRadius)
            maxExclusionRadius = exclusionRadius[iType];
      }

      // the minimum cell size is twice the maxExclusionRadius,
      // but to save memory, we take 2 times that value 
      CellList cellList;
      cellList.setAtomCapacity(system.simulation().atomCapacity());
      cellList.setup(boundary, 2.0*2.0*maxExclusionRadius);

      if (nMolecule > capacity())
         UTIL_THROW("nMolecule > Species.capacity()!"); 

      Simulation& sim = system.simulation();
      for (iMol = 0; iMol < nMolecule; ++iMol) {
         // Add a new molecule to the system
         Molecule &newMolecule= sim.getMolecule(id());
         system.addMolecule(newMolecule);

         // Try placing atoms
         bool moleculeHasBeenPlaced = false;
         for (int iAttempt = 0; iAttempt< maxPlacementAttempts_; iAttempt++) {
            // Place first atom
            Vector pos;
            system.boundary().randomPosition(system.simulation().random(),pos);
            Atom &thisAtom = newMolecule.atom(0);
 
            // check if the first atom can be placed at the new position
            CellList::NeighborArray neighbors;
            cellList.getNeighbors(pos, neighbors);
            int nNeighbor = neighbors.size();
            bool canBePlaced = true;
            for (int j = 0; j < nNeighbor; ++j) {
               Atom *jAtomPtr = neighbors[j];
         
               double r = sqrt(system.boundary().distanceSq(
                                        jAtomPtr->position(), pos));
               if (r < (exclusionRadius[thisAtom.typeId()] +
                  exclusionRadius[jAtomPtr->typeId()])) {
                  canBePlaced = false;
                  break;
               }
            } 
            if (canBePlaced)  {
               thisAtom.position() = pos;
               cellList.addAtom(thisAtom);

               // Try to recursively place other atoms
               if (tryPlaceAtom(newMolecule, 0, exclusionRadius, system,
                  cellList, bondPotentialPtr, system.boundary())) {
                  moleculeHasBeenPlaced = true;
                  break;
              } else {
                 cellList.deleteAtom(thisAtom); 
              }
            }
         }
         if (! moleculeHasBeenPlaced) {
           std::ostringstream oss;
           oss <<  "Failed to place molecule " << newMolecule.id();
           UTIL_THROW(oss.str().c_str());
         }

      }

      #if 0
      // Check
      for (int iMol =0; iMol < nMolecule; ++iMol) {
         Molecule::AtomIterator atomIter;
         system.molecule(id(),iMol).begin(atomIter);
         for (; atomIter.notEnd(); ++atomIter) {
            for (int jMol =0; jMol < nMolecule; ++jMol) {
               Molecule::AtomIterator atomIter2;
               system.molecule(id(),jMol).begin(atomIter2);
               for (; atomIter2.notEnd(); ++atomIter2 ) {
                  if (atomIter2->id() != atomIter->id()) {
                     double r = sqrt(boundary.distanceSq(
                        atomIter->position(),atomIter2->position()));
                     if (r < (exclusionRadius[atomIter->typeId()]+
                        exclusionRadius[atomIter2->typeId()])) {
                        std::cout << r << std::endl;
                        UTIL_THROW("ERROR");
                     }
                  }
               }
            }
         }
      }
      #endif

   } 
   
}
