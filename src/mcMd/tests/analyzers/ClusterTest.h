#ifndef MCMD_CLUSTER_TEST_H
#define MCMD_CLUSTER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/analyzers/system/Cluster.h>
#include <mcMd/analyzers/system/ClusterLink.h>
#include <util/containers/DArray.h>

#include <fstream>

using namespace Util;
using namespace McMd;

class ClusterTest : public UnitTest 
{

   DArray<McMd::ClusterLink> links_;
   McMd::Cluster cluster_;

public:

   void setUp()
   {
      int capacity = 10;
      links_.allocate(capacity); 
      for (int i = 0; i < capacity; ++i) {
         links_[i].clear();
      }
   }

   void tearDown()
   {}

   void addLinks() {
      printMethod(TEST_FUNC);
      Cluster v;
      cluster_.clear();
      cluster_.setId(35);
      TEST_ASSERT(links_[6].clusterId() == -1);
      cluster_.addLink(links_[6]);
      TEST_ASSERT(links_[6].clusterId() == 35);
      TEST_ASSERT(links_[6].next() == 0);
      TEST_ASSERT( cluster_.size() == 1);
      cluster_.addLink(links_[3]);
      TEST_ASSERT(links_[3].clusterId() == 35);
      TEST_ASSERT(links_[3].next() == &links_[6]);
      TEST_ASSERT( cluster_.size() == 2);
      cluster_.addLink(links_[8]);
      TEST_ASSERT(links_[8].clusterId() == 35);
      TEST_ASSERT(links_[8].next() == &links_[3]);
      TEST_ASSERT(cluster_.size() == 3);
      TEST_ASSERT(cluster_.isValid());
      cluster_.setId(48);
      TEST_ASSERT(links_[6].clusterId() == 48);
      TEST_ASSERT(links_[3].clusterId() == 48);
      TEST_ASSERT(links_[8].clusterId() == 48);
      TEST_ASSERT(cluster_.size() == 3);
      TEST_ASSERT(cluster_.isValid());
      cluster_.clear();
      TEST_ASSERT(cluster_.id() == -1);
      TEST_ASSERT(cluster_.head() == 0);
      TEST_ASSERT(cluster_.size() == 0);
      TEST_ASSERT(cluster_.isValid());
      for (int i = 0; i < links_.capacity(); ++i) {
         TEST_ASSERT(links_[i].clusterId() == -1);
         TEST_ASSERT(links_[i].next() == 0);
      }
   }

};

TEST_BEGIN(ClusterTest)
TEST_ADD(ClusterTest, addLinks)
TEST_END(ClusterTest)


#endif
