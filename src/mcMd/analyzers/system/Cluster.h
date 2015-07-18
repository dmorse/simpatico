#ifndef MCMD_CLUSTER_H
#define MCMD_CLUSTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ClusterLink.h"

namespace McMd
{

   /**
   * Cluster of molecules.
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
      * Get the number of molecules/links in the cluster.
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
