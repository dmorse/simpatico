#ifndef LINK_MASTER_TEST_H
#define LINK_MASTER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/links/LinkMaster.h>
#include <mcMd/links/LinkEvents.h>
#include <util/containers/RArray.h>
#include <mcMd/chemistry/Atom.h>
#include <util/util/Observer.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

class LinkObserver : public Observer<LinkAddEvent>, public Observer<LinkRemoveEvent>
{

public:

  LinkObserver()
   : nAdd_(0),
     nRemove_(0)
  {}

  ~LinkObserver()
  {}

  void initialize(LinkMaster& linkMaster)
  {
     linkMaster.Notifier<LinkAddEvent>::registerObserver(*this);
     linkMaster.Notifier<LinkRemoveEvent>::registerObserver(*this);
  }

  virtual void update(const LinkAddEvent& event)
  { 
     ++nAdd_;
     //std::cout << "Added" << std::endl; 
  }

  virtual void update(const LinkRemoveEvent& event)
  { 
    ++nRemove_;
    //std::cout << "Removed" << std::endl; 
  }

  int nAdd() 
  { return nAdd_; }

  int nRemove() 
  { return nRemove_; }

private:

  int nAdd_ ;
  int nRemove_ ;

};

class LinkMasterTest : public UnitTest 
{

public:

   void setUp()
   {
      atomCapacity_ = 100;
      Atom::allocate(atomCapacity_, atoms_);

      std::ifstream file("in/LinkMaster");
      linkMaster_.readParam(file);
      file.close();

      observer_.initialize(linkMaster_);

   }

   void tearDown()
   {
      Atom::deallocate();
   }

   void testReadParam() 
   {
      printMethod(TEST_FUNC);

      std::cout << std::endl;
      linkMaster_.writeParam(std::cout);

      TEST_ASSERT(linkMaster_.nLink() == 0);
      linkMaster_.isValid();
   }

   void testAdd() 
   {
      printMethod(TEST_FUNC);
      const Link* linkPtr0;
      const Link* linkPtr1;
      const LinkMaster::AtomLinkSet* setPtr0; 
      const LinkMaster::AtomLinkSet* setPtr1; 

      TEST_ASSERT(linkMaster_.nLink() == 0);
      linkMaster_.addLink(atoms_[5], atoms_[21], 0);
      TEST_ASSERT(linkMaster_.nLink() == 1);
      setPtr0  = &linkMaster_.atomLinkSet(atoms_[5]);
      setPtr1  = &linkMaster_.atomLinkSet(atoms_[21]);
      TEST_ASSERT(setPtr0->size() == 1);
      TEST_ASSERT(setPtr1->size() == 1);
      linkPtr0 = &linkMaster_.link(0);
      TEST_ASSERT(&linkPtr0->atom0() == &atoms_[5]);
      TEST_ASSERT(&linkPtr0->atom1() == &atoms_[21]);
      TEST_ASSERT(setPtr0->isElement(*linkPtr0));
      TEST_ASSERT(setPtr1->isElement(*linkPtr0));
      linkMaster_.isValid();

      linkMaster_.addLink(atoms_[9], atoms_[5], 0);
      TEST_ASSERT(linkMaster_.nLink() == 2);
      setPtr0 = &linkMaster_.atomLinkSet(atoms_[9]);
      setPtr1 = &linkMaster_.atomLinkSet(atoms_[5]);
      TEST_ASSERT(setPtr0->size() == 1);
      TEST_ASSERT(setPtr1->size() == 2);
      linkPtr0 = &linkMaster_.link(0);
      linkPtr1 = &linkMaster_.link(1);
      TEST_ASSERT(&linkPtr0->atom0() == &atoms_[5]);
      TEST_ASSERT(&linkPtr0->atom1() == &atoms_[21]);
      TEST_ASSERT(&linkPtr1->atom0() == &atoms_[9]);
      TEST_ASSERT(&linkPtr1->atom1() == &atoms_[5]);
      TEST_ASSERT(setPtr0->isElement(*linkPtr1));
      TEST_ASSERT(setPtr1->isElement(*linkPtr1));
      TEST_ASSERT(setPtr1->isElement(*linkPtr0));
      TEST_ASSERT(&(*setPtr0)[0] == linkPtr1);
      TEST_ASSERT(&(*setPtr1)[0] == linkPtr0);
      TEST_ASSERT(&(*setPtr1)[1] == linkPtr1);
      linkMaster_.isValid();

      TEST_ASSERT(observer_.nAdd() == 2);
      TEST_ASSERT(observer_.nRemove() == 0);

   }

   void testAddRemove1() 
   {
      printMethod(TEST_FUNC);
      const Link* linkPtr0;
      const Link* linkPtr1;
      const LinkMaster::AtomLinkSet* setPtr0; 
      const LinkMaster::AtomLinkSet* setPtr1; 
      const LinkMaster::AtomLinkSet* setPtr2; 

      TEST_ASSERT(linkMaster_.nLink() == 0);

      linkMaster_.addLink(atoms_[5], atoms_[21], 0);
      linkMaster_.isValid();
      linkMaster_.addLink(atoms_[9], atoms_[5], 0);
      linkMaster_.isValid();
      linkMaster_.addLink(atoms_[21], atoms_[9], 0);
      linkMaster_.isValid();
      linkMaster_.removeLink(1);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 2);

      setPtr0 = &linkMaster_.atomLinkSet(atoms_[5]);
      setPtr1 = &linkMaster_.atomLinkSet(atoms_[9]);
      setPtr2 = &linkMaster_.atomLinkSet(atoms_[21]);
      TEST_ASSERT(setPtr0->size() == 1);
      TEST_ASSERT(setPtr1->size() == 1);
      TEST_ASSERT(setPtr2->size() == 2);
      linkPtr0 = &linkMaster_.link(0);
      linkPtr1 = &linkMaster_.link(1);
      TEST_ASSERT(&linkPtr0->atom0() == &atoms_[5]);
      TEST_ASSERT(&linkPtr0->atom1() == &atoms_[21]);
      TEST_ASSERT(&linkPtr1->atom0() == &atoms_[21]);
      TEST_ASSERT(&linkPtr1->atom1() == &atoms_[9]);
      TEST_ASSERT(setPtr0->isElement(*linkPtr0));
      TEST_ASSERT(setPtr1->isElement(*linkPtr1));
      TEST_ASSERT(setPtr2->isElement(*linkPtr0));
      TEST_ASSERT(setPtr2->isElement(*linkPtr1));
      TEST_ASSERT(&(*setPtr0)[0] == linkPtr0);
      TEST_ASSERT(&(*setPtr1)[0] == linkPtr1);
      TEST_ASSERT(&(*setPtr2)[0] == linkPtr0);
      TEST_ASSERT(&(*setPtr2)[1] == linkPtr1);

      TEST_ASSERT(observer_.nAdd() == 3);
      TEST_ASSERT(observer_.nRemove() == 1);

   }

   void testAddRemove2() 
   {
      printMethod(TEST_FUNC);
      const Link* linkPtr0;
      const Link* linkPtr1;
      const LinkMaster::AtomLinkSet* setPtr0; 
      const LinkMaster::AtomLinkSet* setPtr1; 
      const LinkMaster::AtomLinkSet* setPtr2; 

      TEST_ASSERT(linkMaster_.nLink() == 0);

      linkMaster_.addLink(atoms_[5], atoms_[21], 0);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 1);

      linkMaster_.addLink(atoms_[9], atoms_[5], 0);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 2);

      linkMaster_.addLink(atoms_[21], atoms_[9], 0);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 3);

      linkMaster_.removeLink(0);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 2);

      setPtr0 = &linkMaster_.atomLinkSet(atoms_[5]);
      setPtr1 = &linkMaster_.atomLinkSet(atoms_[9]);
      setPtr2 = &linkMaster_.atomLinkSet(atoms_[21]);
      TEST_ASSERT(setPtr0->size() == 1);
      TEST_ASSERT(setPtr1->size() == 2);
      TEST_ASSERT(setPtr2->size() == 1);
      linkPtr0 = &linkMaster_.link(0);
      linkPtr1 = &linkMaster_.link(1);
      TEST_ASSERT(&linkPtr0->atom0() == &atoms_[21]);
      TEST_ASSERT(&linkPtr0->atom1() == &atoms_[9]);
      TEST_ASSERT(&linkPtr1->atom0() == &atoms_[9]);
      TEST_ASSERT(&linkPtr1->atom1() == &atoms_[5]);
      TEST_ASSERT(setPtr1->isElement(*linkPtr0));
      TEST_ASSERT(setPtr2->isElement(*linkPtr0));
      TEST_ASSERT(setPtr0->isElement(*linkPtr1));
      TEST_ASSERT(setPtr1->isElement(*linkPtr1));
      TEST_ASSERT(&(*setPtr0)[0] == linkPtr1);
      TEST_ASSERT(&(*setPtr1)[0] == linkPtr1);
      TEST_ASSERT(&(*setPtr1)[1] == linkPtr0);
      TEST_ASSERT(&(*setPtr2)[0] == linkPtr0);

      linkMaster_.removeLink(1);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 1);

      TEST_ASSERT(setPtr0->size() == 0);
      TEST_ASSERT(setPtr1->size() == 1);
      TEST_ASSERT(setPtr2->size() == 1);
      linkPtr0 = &linkMaster_.link(0);
      TEST_ASSERT(&linkPtr0->atom0() == &atoms_[21]);
      TEST_ASSERT(&linkPtr0->atom1() == &atoms_[9]);
      TEST_ASSERT(setPtr1->isElement(*linkPtr0));
      TEST_ASSERT(setPtr2->isElement(*linkPtr0));
      TEST_ASSERT(&(*setPtr1)[0] == linkPtr0);
      TEST_ASSERT(&(*setPtr2)[0] == linkPtr0);

      TEST_ASSERT(observer_.nAdd() == 3);
      TEST_ASSERT(observer_.nRemove() == 2);

   }

   void testIterator() 
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(linkMaster_.nLink() == 0);

      linkMaster_.addLink(atoms_[5], atoms_[21], 0);
      linkMaster_.isValid();
      linkMaster_.addLink(atoms_[9], atoms_[5], 0);
      linkMaster_.isValid();
      linkMaster_.addLink(atoms_[21], atoms_[9], 0);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 3);

      LinkMaster::LinkIterator iterator;
      linkMaster_.begin(iterator);
      int i = 0;
      for ( ; !iterator.atEnd(); ++iterator) {
          i++;
      }
      TEST_ASSERT(i == 3);

   }

   #if 0
   void testRandom() 
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(linkMaster_.nLink() == 0);

      linkMaster_.addLink(atoms_[5], atoms_[21], 0);
      linkMaster_.isValid();
      linkMaster_.addLink(atoms_[9], atoms_[5], 0);
      linkMaster_.isValid();
      linkMaster_.addLink(atoms_[21], atoms_[9], 0);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 3);

      Link* link = linkMaster_.

   }
   #endif

   void testClear() 
   {

      printMethod(TEST_FUNC);
      const LinkMaster::AtomLinkSet* setPtr0; 
      const LinkMaster::AtomLinkSet* setPtr1; 
      const LinkMaster::AtomLinkSet* setPtr2; 

      TEST_ASSERT(linkMaster_.nLink() == 0);

      linkMaster_.addLink(atoms_[5], atoms_[21], 0);
      linkMaster_.addLink(atoms_[9], atoms_[5], 0);
      linkMaster_.addLink(atoms_[21], atoms_[9], 0);
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 3);

      linkMaster_.clear();
      linkMaster_.isValid();
      TEST_ASSERT(linkMaster_.nLink() == 0);

      setPtr0 = &linkMaster_.atomLinkSet(atoms_[5]);
      setPtr1 = &linkMaster_.atomLinkSet(atoms_[9]);
      setPtr2 = &linkMaster_.atomLinkSet(atoms_[21]);
      TEST_ASSERT(setPtr0->size() == 0);
      TEST_ASSERT(setPtr1->size() == 0);
      TEST_ASSERT(setPtr2->size() == 0);

   }

private:

   LinkMaster   linkMaster_;
   
   LinkObserver observer_;

   RArray<Atom> atoms_;

   int atomCapacity_;

};

TEST_BEGIN(LinkMasterTest)
TEST_ADD(LinkMasterTest, testReadParam)
TEST_ADD(LinkMasterTest, testAdd)
TEST_ADD(LinkMasterTest, testAddRemove1)
TEST_ADD(LinkMasterTest, testAddRemove2)
TEST_ADD(LinkMasterTest, testIterator)
TEST_ADD(LinkMasterTest, testClear)
TEST_END(LinkMasterTest)

#endif
