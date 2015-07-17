#ifndef MCMD_CLUSTER_IDENTIFIER_CPP
#define MCMD_CLUSTER_IDENTIFIER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ClusterIdentifier.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   ClusterIdentifier::ClusterIdentifier(System& system)
    : molecules_(),
      clusters_(),
      workStack_(),
      cellList_(),
      systemPtr_(&system),
      speciesPtr_(),
      speciesId_(),
      atomTypeId_(),
      cutoff_()
   {}

   /*
   * Initial setup.
   */
   void
   ClusterIdentifier::setup(int speciesId, int atomTypeId, double cutoff)
   {
      speciesId_ = speciesId;
      atomTypeId_ = atomTypeId;
      cutoff_ = cutoff;
      speciesPtr_ = &system().simulation().species(speciesId);
      int nMolecule = speciesPtr_->capacity();
      molecules_.allocate(nMolecule);
      clusters_.reserve(64);
      int nAtom = nMolecule * speciesPtr_->nAtom();
      cellList_.allocate(nAtom, system().boundary(), cutoff_);
   }

   /*
   * Process one molecule.
   */
   void ClusterIdentifier::processNextMolecule(Cluster& cluster)
   {
      ClusterMolecule* thisMolPtr = &workStack_.pop();
      Molecule::AtomIterator atomIter;
      CellList::NeighborArray atomPtrArray;
      ClusterMolecule* otherMolPtr;
      int thisMolId = thisMolPtr->self().id();
      int otherMolId;
      for (thisMolPtr->self().begin(atomIter); atomIter.notEnd(); ++atomIter) {
         if (atomIter->typeId() == atomTypeId_) {
            cellList_.getNeighbors(atomIter->position(), atomPtrArray);
            for (int i = 0; i < atomPtrArray.size(); i++) {
               otherMolId = atomPtrArray[i]->molecule().id();
               if (otherMolId != thisMolId) {
                  otherMolPtr = &molecules_[otherMolId];
                  if (otherMolPtr->clusterId() == -1) {
                     cluster.addMolecule(*otherMolPtr);
                     workStack_.push(*otherMolPtr);
                  } else
                  if (otherMolPtr->clusterId() != cluster.id()) {
                     UTIL_THROW("Cluster Clash!");
                  }
               }
            }
         }
      }
   }

   /*
   * Identify all clusters in the system.
   */
   void ClusterIdentifier::identifyClusters()
   {
      // Clear molecules and clusters array
      for (int i = 0; i < molecules_.capacity(); ++i) {
          molecules_[i].clear();
      }
      clusters_.clear();

      // Build a cellList of relevant atoms
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      cellList_.clear();
      system().begin(speciesId_, molIter);
      for ( ; molIter.notEnd(); ++molIter) {
         for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            if (atomIter->typeId() == atomTypeId_) {
               system().boundary().shift(atomIter->position());
               cellList_.addAtom(*atomIter);
            }
         }
      }

      // Identify clusters
      Cluster* clusterPtr;
      ClusterMolecule* molPtr;
      int clusterId = 0;
      system().begin(speciesId_, molIter);
      for ( ; molIter.notEnd(); ++molIter) {
         molPtr = &(molecules_[molIter->id()]);
         // molPtr->clear();
         molPtr->self_ = molIter.get();
         if (molPtr->clusterId() == -1) {
            clusters_.resize(clusterId+1);
            clusterPtr = &clusters_[clusterId];
            clusterPtr->clear();
            clusterPtr->setId(clusterId);
            clusterPtr->addMolecule(*molPtr);
            workStack_.push(*molPtr);
            while (workStack_.size() > 0) {
               processNextMolecule(*clusterPtr);
            }
            clusterId++;
         }
      }

   }

}
#endif
