#ifndef MCMD_CLUSTER_IDENTIFIER_H
#define MCMD_CLUSTER_IDENTIFIER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>       // base class template
#include <mcMd/simulation/System.h>              // class template parameter
#include <mcMd/neighbor/CellList.h>              // member
#include <util/containers/DSArray.h>             // member template
#include <util/containers/GArray.h>              // member template
#include <util/containers/GStack.h>              // member template

#include <cstdio>
#include <cstring> 

namespace McMd
{
   using namespace Util;

   class Species;

   /**
   * This class is intended to identify Clusters in polymeric systems.
   */
   class ClusterIdentifier 
   {
   
   public:

      /**
      * Molecule in a cluster.
      */
      struct ClusterMolecule
      {
         Molecule* self;
         ClusterMolecule* next;
         int clusterId;
      };

      /**
      * Cluster.
      */
      struct Cluster
      {
         int id;
         int size;
         ClusterMolecule* head;
      };

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      ClusterIdentifier(System &system);
   
      /** 
      * Clear accumulator.
      */
      virtual void setup(int speciesId, int coreId, double cutoff);
   
      /**
      * Find all clusters.
      */
      void identifyClusters();
 
      int nMolecule() const
      {  return molecules_.size(); }

      int nCluster() const
      {  return clusters_.size(); }

      ClusterMolecule& molecule(int i)
      {  return molecules_[i]; }

      Cluster& cluster(int i)
      {  return clusters_[i]; }

private:

      /// Array of cluster molecule objects.
      DSArray<ClusterMolecule>  molecules_;

      /// Array of length of different Clusters
      GArray<Cluster> clusters_;

      /// Work stack.
      GStack<ClusterMolecule>  workStack_;

      /// CellList of specified atoms of the species of interest
      CellList cellList_;

      /// Pointer to relevant species
      System* systemPtr_;

      /// Pointer to relevant species
      Species* speciesPtr_;

      /// Molecule species type id
      int speciesId_;

      /// Atom typeId of core atoms
      int atomTypeId_;

      /// Cutoff distance for touching cores
      double cutoff_;

      System& system()
      {  return *systemPtr_; }

      /*
      * Process top molecule in the workStack.
      */
      void processNextMolecule(int clusterId);

   };
}
#endif
