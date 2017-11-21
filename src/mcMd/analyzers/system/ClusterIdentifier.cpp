#ifndef MCMD_CLUSTER_IDENTIFIER_CPP
#define MCMD_CLUSTER_IDENTIFIER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ClusterIdentifier.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>
#include <simp/boundary/Boundary.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

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
   ClusterIdentifier::initialize(int speciesId, int atomTypeId, double cutoff)
   {
      // Set member variables
      speciesId_ = speciesId;
      atomTypeId_ = atomTypeId;
      cutoff_ = cutoff;

      // Allocate memory
      Species* speciesPtr = &system().simulation().species(speciesId);
      int moleculeCapacity = speciesPtr->capacity();
      links_.allocate(moleculeCapacity);
      clusters_.reserve(64);
      int atomCapacity = system().simulation().atomCapacity();
      cellList_.setAtomCapacity(atomCapacity);

      // Note: We must set the cellist atom capacity to the total atom capacity,
      // even though we are only interested in clusters of one species, because 
      // the celllist atom capacity sets the maximum allowed atom index value.
   }

   /*
   * Pop and process the next molecule on the workStack (private).
   */
   void ClusterIdentifier::processNextMolecule(Cluster& cluster)
   {
      CellList::NeighborArray neighborArray;
      Molecule::AtomIterator atomIter;
      ClusterLink* thisLinkPtr;
      ClusterLink* otherLinkPtr;
      Atom* otherAtomPtr;
      Boundary& boundary = system().boundary();
      double cutoffSq = cutoff_*cutoff_;
      double rsq;
      int thisMolId, otherMolId, otherClusterId;

      // Pop this molecule off the stack
      thisLinkPtr = &workStack_.pop();
      thisMolId = thisLinkPtr->molecule().id();
      if (thisLinkPtr->clusterId() != cluster.id()) {
         UTIL_THROW("Top ClusterLink not marked with this cluster id");
      }

      /*
      * Loop over atoms of this molecule.
      * For each atom of type atomTypeId_, find neighboring molecules.
      * Add each new neighbor to the cluster, and to the workStack.
      */
      thisLinkPtr->molecule().begin(atomIter); 
      for ( ; atomIter.notEnd(); ++atomIter) {
         if (atomIter->typeId() == atomTypeId_) {
            cellList_.getNeighbors(atomIter->position(), neighborArray);
            for (int i = 0; i < neighborArray.size(); i++) {
               otherAtomPtr = neighborArray[i];
               otherMolId = otherAtomPtr->molecule().id();
               if (otherMolId != thisMolId) {
                  rsq = boundary.distanceSq(atomIter->position(),
                                            otherAtomPtr->position());
                  if (rsq < cutoffSq) {
                     otherLinkPtr = &(links_[otherMolId]);
                     assert(&otherLinkPtr->molecule() ==
                            &otherAtomPtr->molecule());
                     otherClusterId = otherLinkPtr->clusterId();
                     if (otherClusterId == -1) {
                        cluster.addLink(*otherLinkPtr);
                        workStack_.push(*otherLinkPtr);
                     } else
                     if (otherClusterId != cluster.id()) {
                        UTIL_THROW("Cluster Clash!");
                     }
                  }
               }
            } // neighbor atoms
         }
      } // atoms

   }

   /*
   * Identify all clusters in the system.
   */
   void ClusterIdentifier::identifyClusters()
   {

      // Initialize all data structures:
      // Setup a grid of empty cells
      cellList_.setup(system().boundary(), cutoff_);
      // Clear clusters array and all links
      clusters_.clear();
      for (int i = 0; i < links_.capacity(); ++i) {
         links_[i].clear();
      }

      // Build the cellList, associate Molecule with ClusterLink.
      // Iterate over molecules of species speciesId_
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      system().begin(speciesId_, molIter);
      for ( ; molIter.notEnd(); ++molIter) {

         // Associate this Molecule with a ClusterLink
         links_[molIter->id()].setMolecule(*molIter.get());

         // Add atoms of type = atomTypeId_ to the CellList
         for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            if (atomIter->typeId() == atomTypeId_) {
               system().boundary().shift(atomIter->position());
               cellList_.addAtom(*atomIter);
            }

         }
      }

      // Identify all clusters
      Cluster* clusterPtr;
      ClusterLink* linkPtr;
      int clusterId = 0;
      system().begin(speciesId_, molIter);
      for ( ; molIter.notEnd(); ++molIter) {

         // Find the link with same index as this molecule
         linkPtr = &(links_[molIter->id()]);
         assert (&(linkPtr->molecule()) == molIter.get());

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

      // Validity check - throws exception on failure.
      isValid();
   }

   bool ClusterIdentifier::isValid() const
   {
      // Check clusters
      int nCluster = clusters_.size();
      int nMolecule = 0;
      for (int i = 0; i < nCluster; ++i) {
         if (!clusters_[i].isValid()) {
            UTIL_THROW("Invalid cluster");
         }
         nMolecule += clusters_[i].size();
      }
      if (nMolecule != systemPtr_->nMolecule(speciesId_)) {
         UTIL_THROW("Error in number of molecules");
      }
   
      // Check molecules and links
      ClusterLink const * linkPtr;
      System::ConstMoleculeIterator molIter;
      systemPtr_->begin(speciesId_, molIter);
      for ( ; molIter.notEnd(); ++molIter) {

         linkPtr = &(links_[molIter->id()]);
         if (&(linkPtr->molecule()) != molIter.get()) {
            UTIL_THROW("Link without correct molecule association");
         }
         if (linkPtr->clusterId() == -1) {
            UTIL_THROW("Unmarked molecule");
         }
      }

      // Normal return (no errors)
      return true;
   }

}
#endif
