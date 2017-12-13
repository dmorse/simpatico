#ifndef MCMD_CLUSTER_H
#define MCMD_CLUSTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ClusterLink.h"
#include <util/space/Tensor.h>
#include <util/space/Vector.h>
#include <util/global.h>
#include <mcMd/simulation/System.h>                   // base class templ param

namespace McMd
{
   class System;
   
   /**
   * Cluster of molecules.
   *
   * A Cluster is implemented as a linked list of ClusterLink
   * objects, each of which is associated with a molecule. A
   * Cluster has a pointer to the first ClusterLink and a size 
   * member that counts the number of links (or molecules), but 
   * does not own the associated ClusterLink or Molecule 
   * objects.
   */
   class Cluster
   {
   
   public:

      /**
      * Constructor.
      */
      Cluster();

      /**
      * Destructor.
      */
      ~Cluster();

      /**
      * Set cluster to empty.
      */
      void clear();

      /**
      * Set cluster identifier.
      */
      void setId(int id);

      /**
      * Add a link to the list.
      *
      * \param link ClusterLink associated with a Molecule.
      */
      void addLink(ClusterLink& link);

      /**
      * Get the cluster id.
      */
      int id() const
      {  return id_; }

      /**
      * Get the number of molecules or links in the cluster.
      */
      int size() const
      {  return size_; }

      /**
      * Get a pointer to the first link in the linked list.
      *
      * Returns 0 pointer if cluster is empty.
      */
      ClusterLink* head() const
      {  return head_; }

      /**
      * Return true if valid, false otherwise.
      */
      bool isValid() const;

      /**
      * Return the cluster COM
      */
      Vector clusterCOM( int atomTypeInCluster, Boundary const & boundary);

      /**
      * Return the cluster radius of gyration tensor
      */
      Tensor momentTensor(int atomTypeInCluster, Boundary const & boundary);

   private:

      // Integer identifier for this cluster.
      int id_;

      // Number of molecules (or links) in this cluster.
      int size_;

      // Pointer to first link in singly-linked list.
      ClusterLink* head_;

   };

}
#endif
