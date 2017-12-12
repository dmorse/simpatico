/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Cluster.h"
#include <util/global.h>

using namespace Util;

namespace McMd
{

   Cluster::Cluster() :
      id_(-1),
      size_(0),
      head_(0)
   {}

   Cluster::~Cluster() 
   {}

   void Cluster::clear()
   {
      ClusterLink* next = 0;
      ClusterLink* ptr = head_;
      while (ptr) {
         next = ptr->next();
         ptr->clear();
         ptr = next;
      }
      id_ = -1;
      size_ = 0;
      head_ = 0;
   }

   void  Cluster::setId(int id)
   { 
      if (id < 0) {
         UTIL_THROW("Cluster id cannot be negative");
      }
      id_ = id; 
      ClusterLink* ptr = head_;
      while (ptr) {
         ptr->clusterId_ = id;
         ptr = ptr->next();
      }
   }

   void Cluster::addLink(ClusterLink& link)
   {
      link.clusterId_ = id_;
      link.next_ = head_;
      head_ = &link;
      ++size_;
   }

   bool Cluster::isValid() const
   {
      int n = 0;
      ClusterLink* ptr = head_;
      while (ptr) {
         ++n;
         if (ptr->clusterId() != id_) {
            return false;
         }
         ptr = ptr->next();
      }
      if (n != size_) {
         return false;
      }
      return true;
   }

   Vector Cluster::clusterCOM(int atomTypeInCluster, Boundary const & boundary)
   {
     ClusterLink* thisClusterStart;
     ClusterLink* next; 
     Vector com;
     Vector dr;
     com.zero();
     thisClusterStart = head();
     Vector r0 = (thisClusterStart->molecule()).atom(0).position();
     Molecule thisMolecule;
     Molecule::ConstAtomIterator atomIter; 
     int nAtomsInCluster = 0;     
       

     com.zero();
     while(thisClusterStart) {
       next = thisClusterStart->next();
       thisMolecule = thisClusterStart->molecule();
       thisMolecule.begin(atomIter);
       for( ; atomIter.notEnd(); ++atomIter) {
         if (atomIter->typeId() == atomTypeInCluster) {
           boundary.distanceSq(atomIter->position(),r0,dr);
           com += dr;
           nAtomsInCluster += 1;
         }
       }
       thisClusterStart = next;
     }
     com /= nAtomsInCluster;
     com += r0;
     boundary.shift(com);
     return com;
   }

   Tensor Cluster::momentTensor(int atomTypeInCluster, Boundary const & boundary)
   {
     Vector com = clusterCOM( atomTypeInCluster, boundary);
     Tensor rgTensor;
     rgTensor.zero();
     ClusterLink* thisClusterStart;
     ClusterLink* next; 
     thisClusterStart = head();
     Molecule thisMolecule;
     Molecule::ConstAtomIterator atomIter;
     Vector dr;
     Tensor rgDyad;
     int nAtomsInCluster = 0;     
     while(thisClusterStart) {
       next = thisClusterStart->next();
       thisMolecule = thisClusterStart->molecule();
       thisMolecule.begin(atomIter);
       for( ; atomIter.notEnd(); ++atomIter) {
         if (atomIter->typeId() == atomTypeInCluster) {
           nAtomsInCluster += 1;
           boundary.distanceSq(atomIter->position(), com,dr);
           rgTensor += rgDyad.dyad(dr,dr);
         }
       }
       thisClusterStart = next;
     }
     rgTensor /= nAtomsInCluster;
     return rgTensor;
   }

}
