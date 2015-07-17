#ifndef MCMD_CLUSTER_MOLECULE_H
#define MCMD_CLUSTER_MOLECULE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Molecule;
   class Cluster;
   class ClusterIdentifier;

   /**
   * Molecule in a cluster.
   */
   struct ClusterMolecule
   {

   public:

      void clear()
      {
         self_ = 0;  
         next_ = 0;  
         clusterId_ = -1;  
      }

      Molecule& self() const
      {  return *self_; }

      const ClusterMolecule& next() const
      {  return *next_; }

      int clusterId() const
      {  return clusterId_; }

   private:

      Molecule* self_;
      ClusterMolecule* next_;
      int clusterId_;

      friend class Cluster;
      friend class ClusterIdentifier;

   };

}
#endif
