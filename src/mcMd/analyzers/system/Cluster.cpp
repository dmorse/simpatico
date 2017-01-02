/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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

}
