#ifndef MCMD_CLUSTER_IDENTIFIER_H
#define MCMD_CLUSTER_IDENTIFIER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>                  // base  class template parameter
#include <mcMd/analyzers/system/Cluster.h>           // member
#include <mcMd/analyzers/system/ClusterLink.h>   // member
#include <mcMd/neighbor/CellList.h>                  // member
#include <util/containers/DArray.h>                  // member template
#include <util/containers/GArray.h>                  // member template
#include <util/containers/GStack.h>                  // member template

#include <cstdio>
#include <cstring> 

namespace McMd
{
   using namespace Util;

   class System;
   class Species;

   /**
   * This class is intended to identify Clusters in polymeric systems.
   */
   class ClusterIdentifier 
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      ClusterIdentifier(System &system);
   
      /**
      * Destructor.
      */
      ~ClusterIdentifier();
   
      /** 
      * Clear accumulator.
      */
      virtual void setup(int speciesId, int coreId, double cutoff);
   
      /**
      * Find all clusters.
      */
      void identifyClusters();
 
      int nCluster() const
      {  return clusters_.size(); }

      ClusterLink& link(int i)
      {  return links_[i]; }

      Cluster& cluster(int i)
      {  return clusters_[i]; }

private:

      /// Array of cluster link objects, indexed by molecule id.
      DArray<ClusterLink>  links_;

      /// Growable array of clusters.
      GArray<Cluster> clusters_;

      /// Work stack of unprocessed cluster links.
      GStack<ClusterLink>  workStack_;

      /// CellList of atoms of the specified species and atom type.
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
      * Process top ClusterLink in the workStack.
      */
      void processNextMolecule(Cluster& cluster);

   };
}
#endif
