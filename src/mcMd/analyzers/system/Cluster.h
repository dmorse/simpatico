#ifndef MCMD_CLUSTER_H
#define MCMD_CLUSTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ClusterMolecule.h"

namespace McMd
{

   /**
   * Cluster of molecules.
   */
   class Cluster
   {
   
   public:

      Cluster();

      void clear();

      void setId(int id);

      void addMolecule(ClusterMolecule& clusterMolecule);

      int id() const
      {  return id_; }

      int size() const
      {  return size_; }

      const ClusterMolecule& head() const
      {  return *head_; }

   private:

      int id_;
      int size_;
      ClusterMolecule* head_;

   };

}
#endif
