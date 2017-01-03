#ifndef MCMD_CLUSTER_MOLECULE_H
#define MCMD_CLUSTER_MOLECULE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Molecule;
   class Cluster;

   /**
   * Molecule in a cluster.
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

      Molecule* moleculePtr_;
      ClusterLink* next_;
      int clusterId_;

      friend class Cluster;

   };

   inline
   void ClusterLink::clear()
   {
      moleculePtr_ = 0;  
      next_ = 0;  
      clusterId_ = -1;  
   }

   inline
   void ClusterLink::setMolecule(Molecule& molecule)
   {  moleculePtr_ = &molecule;  }

   inline
   Molecule& ClusterLink::molecule() const
   {  return *moleculePtr_; }

   inline
   ClusterLink* ClusterLink::next() const
   {  return next_; }

   inline
   int ClusterLink::clusterId() const
   {  return clusterId_; }

}
#endif
