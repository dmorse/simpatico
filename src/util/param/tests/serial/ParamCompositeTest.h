#ifndef PARAM_COMPOSITE_TEST_H
#define PARAM_COMPOSITE_TEST_H

#include <util/global.h>

#include <util/param/ParamComposite.h>
#include <util/param/Factory.h>
#include <util/param/Manager.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>

#include <util/archives/MemoryCounter.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>
#include <fstream>

using namespace Util;

#include "../ParamTestClasses.h"

class ParamCompositeTest : public ParamFileTest<ParamComposite>
{

public:

   void setUp()
   {}

   void testConstructor() 
   {}

   void testAddWrite() 
   {
      printMethod(TEST_FUNC);

      int     value0 = 4;
      long    value1 = 25;
      double  value2 = 3.0;
      int     value3[3];
      value3[0] = 4;
      value3[1] = 5;
      value3[2] = 6;
      double value4[3];
      value4[0] = 4.0;
      value4[1] = 5.0;
      value4[2] = 6.0;
      double value5[2][2];
      value5[0][0] = 4.0;
      value5[0][1] = 5.0;
      value5[1][0] = 6.0;
      value5[1][1] = 7.0;
      DArray<double> value6;
      value6.allocate(4);
      value6[0] = 14.5;
      value6[1] = 15.4;
      value6[2] = 16.3;
      value6[3] = 16.2;
      FArray<double, 4> value7;
      value7[0] = 13.5;
      value7[1] = 14.4;
      value7[2] = 15.3;
      value7[3] = 15.2;

      object().addBegin("ClassName");
      object().add<int>("value0", value0);
      object().add<long>("value1", value1);
      object().add<double>("value2", value2);
      object().addCArray<int>("value3", value3, 3);
      object().addCArray<double>("value4", value4, 3);
      object().addCArray2D<double>("value5", &value5[0][0], 2, 2);
      object().addDArray<double>("value6", value6, 4);
      object().addFArray<double, 4>("value7", value7);
      object().addEnd();

      printEndl();
      object().writeParam(std::cout);
   }

   void testReadWrite() 
   {
      printMethod(TEST_FUNC);
      int     value0;
      long    value1;
      double  value2;
      int     value3[3];
      double  value4[3];
      double  value5[2][2];
      DArray<double> value6;
      value6.allocate(4);
      Vector    value7;
      IntVector value8;
      DMatrix<double>  value9;
      value9.allocate(2, 2);
      AManager  manager;

      openFile("in/ParamComposite");

      //object().setEcho();
      object().readBegin(file(), "ClassName");
      object().read<int>(file(), "value0", value0);
      object().read<long>(file(), "value1", value1);
      object().read<double>(file(), "value2", value2);
      object().readBlank(file());
      object().readCArray<int>(file(), "value3", value3, 3);
      object().readCArray<double>(file(), "value4", value4, 3);
      object().readCArray2D<double>(file(), "value5", &value5[0][0], 2, 2);
      object().readDArray<double>(file(), "value6", value6, 4);
      object().read<Vector>(file(), "value7", value7);
      object().read<IntVector>(file(), "value8", value8);
      object().readDMatrix<double>(file(), "value9", value9, 2, 2);
      object().readParamComposite(file(), manager);
      object().readEnd(file());

      printEndl();
      object().writeParam(std::cout);
   }

   void testSerialize() 
   {
      printMethod(TEST_FUNC);

      openFile("in/ParamComposite");

      AComposite original;
      original.readParam(file());

      MemoryCounter  car;
      car << original;
      int size = car.size();

      MemoryOArchive oar;
      oar.allocate(size);
      oar << original;

      MemoryIArchive iar;
      iar = oar;

      AComposite clone;
      iar >> clone;

      printEndl();
      clone.writeParam(std::cout);
   }

};

TEST_BEGIN(ParamCompositeTest)
TEST_ADD(ParamCompositeTest, testConstructor)
TEST_ADD(ParamCompositeTest, testAddWrite)
TEST_ADD(ParamCompositeTest, testReadWrite)
TEST_ADD(ParamCompositeTest, testSerialize)
TEST_END(ParamCompositeTest)

#endif // ifndef PARAM_COMPOSITE_TEST_H
