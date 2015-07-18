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

   /*
   * Constructor.
   */
   ClusterIdentifier::ClusterIdentifier(System& system)
    : links_(),
      clusters_(),
      workStack_(),
      cellList_(),
      systemPtr_(&system),
      speciesId_(),
      atomTypeId_(),
      cutoff_()
   {}

   /*
   * Destructor.
   */
   ClusterIdentifier::~ClusterIdentifier()
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
      Species* speciesPtr = &system().simulation().species(speciesId);
      int nMolecule = speciesPtr->capacity();
      links_.allocate(nMolecule);
      clusters_.reserve(64);
      int nAtom = nMolecule * speciesPtr->nAtom();
      cellList_.allocate(nAtom, system().boundary(), cutoff_);
   }

   /*
   * Process the next molecule on the workStack.
   */
   void ClusterIdentifier::processNextMolecule(Cluster& cluster)
   {
      Molecule::AtomIterator atomIter;
      CellList::NeighborArray neighborArray;
      ClusterLink* thisLinkPtr = &workStack_.pop();
      ClusterLink* otherLinkPtr;
      int thisMolId = thisLinkPtr->molecule().id();
      int otherMolId;
      thisLinkPtr->molecule().begin(atomIter); 
      for ( ; atomIter.notEnd(); ++atomIter) {
         if (atomIter->typeId() == atomTypeId_) {
            cellList_.getNeighbors(atomIter->position(), neighborArray);
            for (int i = 0; i < neighborArray.size(); i++) {
               otherMolId = neighborArray[i]->molecule().id();
               if (otherMolId != thisMolId) {
                  otherLinkPtr = &links_[otherMolId];
                  if (otherLinkPtr->clusterId() == -1) {
                     cluster.addLink(*otherLinkPtr);
                     workStack_.push(*otherLinkPtr);
                  } else
                  if (otherLinkPtr->clusterId() != cluster.id()) {
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

      // Clear clusters array and all links
      clusters_.clear();
      for (int i = 0; i < links_.capacity(); ++i) {
          links_[i].clear();
      }

      // Build a cellList of atoms, and set molecules for links.
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      cellList_.clear();
      system().begin(speciesId_, molIter);
      for ( ; molIter.notEnd(); ++molIter) {

         // Set associated link to point to this molecule
         links_[molIter->id()].setMolecule(*molIter.get());

         // Add atoms of appropriate type to the CellList
         for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            if (atomIter->typeId() == atomTypeId_) {
               system().boundary().shift(atomIter->position());
               cellList_.addAtom(*atomIter);
            }

         }
      }

      // Identify clusters
      Cluster* clusterPtr;
      ClusterLink* linkPtr;
      int clusterId = 0;
      system().begin(speciesId_, molIter);
      for ( ; molIter.notEnd(); ++molIter) {

         // Find link with same index as this molecule
         linkPtr = &(links_[molIter->id()]);
         assert (&(linkPtr->molecule()) = molIter.get());

         // If this link is not in a cluster, begin a new cluster
         if (linkPtr->clusterId() == -1) {

            // Add a new empty cluster to clusters_ array
            clusters_.resize(clusterId+1);
            clusterPtr = &clusters_[clusterId];
            clusterPtr->clear();
            clusterPtr->setId(clusterId);

            // Identify molecules in this cluster
            clusterPtr->addLink(*linkPtr);
            workStack_.push(*linkPtr);
            while (workStack_.size() > 0) {
               processNextMolecule(*clusterPtr);
            }

            clusterId++;
         }

      }

   }

}
#endif
