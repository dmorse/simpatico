#ifndef MPI_PARAM_COMPOSITE_TEST_H
#define MPI_PARAM_COMPOSITE_TEST_H

#include <util/global.h>

#ifdef UTIL_MPI

#include <util/param/ParamComposite.h>
#include <util/param/Factory.h>
#include <util/param/Manager.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include "FileMaster.cpp"

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

#include "../ParamTestClasses.h"

class MpiParamCompositeTest : public ParamFileTest
{

   AComposite object;

public:

   void setUp()
   {}

   void testReadWrite1() 
   {
      printMethod(TEST_FUNC);
      int     value0;
      long    value1;
      double  value2;
      std::string  str;
      int     value3[3];
      double  value4[3];
      double  value5[2][2];
      DArray<double> value6;
      value6.allocate(4);
      Vector    value7;
      IntVector value8;
      DMatrix<double> value9;
      value9.allocate(2,2);
      E         e;
      AManager  manager;

      object.setIoCommunicator(communicator()); 
      //ParamComponent::setEcho(true);

      openFile("in/ParamComposite");

      object.readBegin(file(), "ClassName");
      object.read<int>(file(), "value0", value0);
      object.read<long>(file(), "value1", value1);
      object.read<double>(file(), "value2", value2);
      object.read<std::string>(file(), "str", str);
      object.readCArray<int>(file(), "value3", value3, 3);
      object.readCArray<double>(file(), "value4", value4, 3);
      object.readCArray2D<double>(file(), "value5", value5[0], 2, 2);
      object.readDArray<double>(file(), "value6", value6, 4);
      object.readBlank(file());
      object.read<Vector>(file(), "value7", value7);
      object.read<IntVector>(file(), "value8", value8);
      object.readDMatrix<double>(file(), "value9", value9, 2, 2);
      object.readParamComposite(file(), e);
      object.readParamComposite(file(), manager);
      object.readEnd(file());

      if (mpiRank() == 1) {
         std::cout << std::endl;
         object.writeParam(std::cout);
      }

   }

   void testReadWrite2() 
   {
      printMethod(TEST_FUNC);

      object.setIoCommunicator(communicator()); 
      //ParamComponent::setEcho(true);

      openFile("in/ParamComposite");
      object.readParam(file());

      if (mpiRank() == 1) {
         std::cout << std::endl;
         object.writeParam(std::cout);
      }

   }

   void testReadWrite3() 
   {
      printMethod(TEST_FUNC);

      FileMaster    fileMaster;
      std::ofstream logFile;

      fileMaster.setDirectoryId(mpiRank());
      fileMaster.openInputFile("ParamComposite", file());
      fileMaster.openOutputFile("log", logFile);
      Log::setFile(logFile);

      object.readParam(file());

      Log::file() << std::endl;
      object.writeParam(Log::file());

   }

};


TEST_BEGIN(MpiParamCompositeTest)
TEST_ADD(MpiParamCompositeTest, testReadWrite1)
TEST_ADD(MpiParamCompositeTest, testReadWrite2)
TEST_ADD(MpiParamCompositeTest, testReadWrite3)
TEST_END(MpiParamCompositeTest)

#endif // ifdef UTIL_MPI
#endif // ifndef MPI_PARAM_COMPOSITE_TEST_H
