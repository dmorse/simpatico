#ifndef MCMD_CLUSTER_LINK_H
#define MCMD_CLUSTER_LINK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Molecule;
   class Cluster;

   /**
   * Molecule in a cluster.
   *
   * A cluster is defined by a linked list of ClusterLink objects.
   */
   struct ClusterLink
   {

   public:

      /**
      * Set to default state.
      */
      void clear();

      /**
      * Set pointer to associated molecule.
      */
      void setMolecule(Molecule& molecule);

      /**
      * Get associated molecule by reference.
      */
      Molecule& molecule() const;

      /**
      * Get pointer to next link in linked list.
      */
      ClusterLink* next() const;

      /**
      * Get cluster identifier.
      */
      int clusterId() const;

   private:

      /// Pointer to the associated molecule.
      Molecule* moleculePtr_;

      /// Pointer to the next link.
      ClusterLink* next_;

      /// Integer id of the associated cluster.
      int clusterId_;

      friend class Cluster;

   };

   // Set ClusterLink to default null state.
   inline 
   void ClusterLink::clear()
   {
      moleculePtr_ = 0;  
      next_ = 0;  
      clusterId_ = -1;  
   }

   // Set pointer to associated Molecule.
   inline
   void ClusterLink::setMolecule(Molecule& molecule)
   {  moleculePtr_ = &molecule;  }

   // Get the associated molecule by reference.
   inline
   Molecule& ClusterLink::molecule() const
   {  return *moleculePtr_; }

   // Get a pointer to the next link.
   inline
   ClusterLink* ClusterLink::next() const
   {  return next_; }

   // Get the id of the associated cluster.
   inline
   int ClusterLink::clusterId() const
   {  return clusterId_; }

}
#endif
