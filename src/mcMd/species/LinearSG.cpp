/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/neighbor/CellList.h>
#include <mcMd/chemistry/Molecule.h>
#include "LinearSG.h"
#ifdef UTIL_MPI
#include <mcMd/simulation/McMd_mpi.h>
#endif

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   LinearSG::LinearSG()
    : Linear(),
      SpeciesMutator(),
      beadTypeIds1_(),
      beadTypeIds2_()
   {
      setMutatorPtr(this);
   }

   /*
   * Destructor.
   */
   LinearSG::~LinearSG()
   {}

   /*
   * Call general Species::readParameters() .
   */
   void LinearSG::readParameters(std::istream& in)
   {
      Species::readParameters(in);
   }

   /*
   * Read atom structure and two sets of atom type ids.
   */
   void LinearSG::readSpeciesParam(std::istream& in)
   {
      read<int>(in,"nAtom", nAtom_);
      nBond_ = nAtom_ - 1;
      #ifdef INTER_ANGLE
      nAngle_ = nAtom_ - 2;
      #endif
      #ifdef INTER_DIHEDRAL
      if (nAtom_ > 3)
         nDihedral_ = nAtom_ - 3;
      else
         nDihedral_ = 0;
      #endif
      buildLinear();

      read<Pair <int> >(in, "typeIds", typeIds_);
      beadTypeIds1_.allocate(nAtom_);
      beadTypeIds2_.allocate(nAtom_);
      readDArray<int>(in, "identities1", beadTypeIds1_, nAtom_);
      readDArray<int>(in, "identities2", beadTypeIds2_, nAtom_);

      read<double>(in, "weightRatio", weightRatio_);

      allocateSpeciesMutator(capacity(), 2);

      // Set statistical weights
      double sum = weightRatio_ + 1.0;
      mutator().setWeight(0, weightRatio_/sum);
      mutator().setWeight(1, 1.0/sum);
   }

   /*
   * Return NullIndex for every atom.
   * Set initial typeIds
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateAtomTypeId(int index) const
   { return 0; }

   /*
   * Return 0 for every bond.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateBondTypeId(int index) const
   { return 0; }

   #ifdef INTER_ANGLE
   /*
   * Return 0 for every angle.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateAngleTypeId(int index) const
   { return 0; }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Return 0 for every dihedral.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateDihedralTypeId(int index) const
   { return 0; }
   #endif

   /*
   * Change the type of a specific molecule.
   */
   void LinearSG::setMoleculeState(Molecule& molecule, int stateId)
   {
      int nAtom  = molecule.nAtom();
      for (int i = 0; i < nAtom; ++i) {
         if (stateId == 0) {
            beadIdentities_ = beadTypeIds1_;
         } else {
            beadIdentities_ = beadTypeIds2_;
         }
         molecule.atom(i).setTypeId(beadIdentities_[i]);
      }
      setMoleculeStateId(molecule, stateId);
   }

   /*
   * Recursive function to try to place an atom.
   *
   bool
   LinearSG::tryPlaceAtom(Molecule& molecule,
                          int atomId,
                          DArray<double> exclusionRadius,
                          System& system,
                          CellList &cellList,
                          BondPotential *bondPotentialPtr,
                          const Boundary &boundary)
   {
      Atom& lastAtom = molecule.atom(atomId);
      Atom& thisAtom = molecule.atom(++atomId);
      Random& random = system.simulation().random();
      int bondTypeId = calculateBondTypeId(lastAtom.indexInMolecule());

      bool hasBeenPlaced = false;

      Vector v, newPos;
      double r;
      Atom* jAtomPtr;
      int beta, nNeighbor;
      bool canBePlaced;
      for (int iAttempt = 0; iAttempt < maxPlacementAttempts_; iAttempt++)
      {
         // Draw a random bond vector
         beta = 1;
         random.unitVector(v);
         v *= bondPotentialPtr->randomBondLength(&random, beta, bondTypeId);

         newPos = lastAtom.position();
         newPos += v;
         // shift into simulation cell
         boundary.shift(newPos);

         // check if the atom can be placed at the new position
         CellList::NeighborArray neighbors;
         cellList.getNeighbors(newPos, neighbors);
         nNeighbor = neighbors.size();
         canBePlaced = true;
         for (int j = 0; j < nNeighbor; ++j) {
            jAtomPtr = neighbors[j];

            r = sqrt(boundary.distanceSq(jAtomPtr->position(), newPos));
            if (r < (exclusionRadius[thisAtom.typeId()] +
                     exclusionRadius[jAtomPtr->typeId()])) {
               canBePlaced = false;
               break;
            }
         }
         if (canBePlaced) {
            // Place the particle
            thisAtom.position() = newPos;

            // Add to cell list
            cellList.addAtom(thisAtom);

            // Are we add the end of the chain?
            if (atomId == molecule.nAtom()-1) {
               return true;
            }

            // Recursion step
            if (!tryPlaceAtom(molecule, atomId, exclusionRadius, system,
                              cellList, bondPotentialPtr, boundary) ) {
               // If next monomer cannot be inserted, delete this monomer
               cellList.deleteAtom(thisAtom);
            } else {
               hasBeenPlaced = true;
               break;
            }
         }
      }

      return hasBeenPlaced;
   }

   *
   * Generate random molecules
   
   void
   LinearSG::generateMolecules(int nMolecule,
                               DArray<double> exclusionRadius,
                               System& system,
                               BondPotential *bondPotentialPtr,
                               const Boundary &boundary)
   {
      int iMol;

      // Compute the maximum exclusion radius
      double maxExclusionRadius = 0.0;
      for (int iType = 0; iType < system.simulation().nAtomType(); iType++) {
         if (exclusionRadius[iType] > maxExclusionRadius) {
            maxExclusionRadius = exclusionRadius[iType];
         }
      }

      // Allocate a cell list.
      // The minimum cell size is twice the maxExclusionRadius,
      // but to save memory, we take 2 times that value
      CellList cellList;
      cellList.allocate(system.simulation().atomCapacity(),
                        boundary, 2.0*2.0*maxExclusionRadius);

      if (nMolecule > capacity()) {
         UTIL_THROW("nMolecule > Species.capacity()!");
      }

      Simulation& sim = system.simulation();
      for (iMol = 0; iMol < nMolecule; ++iMol) {
         // Add a new molecule to the system
         Molecule &newMolecule= sim.getMolecule(id());
         system.addMolecule(newMolecule);

         // Try placing atoms
         bool moleculeHasBeenPlaced = false;
         for (int iAttempt = 0; iAttempt < maxPlacementAttempts_; iAttempt++) {
            // Place first atom
            Vector pos;
            system.boundary().randomPosition(system.simulation().random(), pos);
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
                  setMoleculeState(newMolecule , 0);
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
*/
}
