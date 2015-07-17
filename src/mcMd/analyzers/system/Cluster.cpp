/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Cluster.h"

namespace McMd
{

   Cluster::Cluster() :
      id_(-1),
      size_(0),
      head_(0)
   {}

   void Cluster::clear()
   {
      id_ = -1;
      size_ = 0;
      head_ = 0;
   }

   void  Cluster::setId(int id)
   {
      id_ = id;
   }

   void Cluster::addMolecule(ClusterMolecule& clusterMolecule)
   {
      clusterMolecule.clusterId_ = id_;
      clusterMolecule.next_ = head_;
      head_ = &clusterMolecule;
      ++size_;
   }

}
