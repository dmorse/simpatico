#ifndef PARAM_TEST_H
#define PARAM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/param/ScalarParam.h>
#include <util/param/CArrayParam.h>
#include <util/param/FArrayParam.h>
#include <util/param/CArray2DParam.h>
#include <util/param/DMatrixParam.h>

using namespace Util;

class ParamTest : public UnitTest 
{

public:

   void setUp()
   { setVerbose(2); }

   void tearDown()
   {}

   void testParamIntConstructor() {
      printMethod(TEST_FUNC);
      int        value = 4;
      Parameter* param;
      param = new ScalarParam<int>("MyLabel", value);
      delete param;
   }

   void testParamIntWrite() {
      printMethod(TEST_FUNC);
      int        value = 4;
      Parameter *param;
      param = new ScalarParam<int>("MyLabel", value);

      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testParamIntRead() {
      printMethod(TEST_FUNC);
      int        value;
      Parameter *param;
      param = new ScalarParam<int>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamInt", in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testParamDoubleWrite() {
      printMethod(TEST_FUNC);
      double value = 4.0;
      Parameter     *param;
      param = new ScalarParam<double>("MyLabel", value);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testParamDoubleRead() {
      printMethod(TEST_FUNC);
      double value;
      Parameter *param;
      param = new ScalarParam<double>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamDouble", in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testParamStringWrite() {
      printMethod(TEST_FUNC);
      std::string value = "stringy";
      Parameter *param;
      param = new ScalarParam<std::string>("MyLabel", value);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testParamStringRead() {
      printMethod(TEST_FUNC);
      std::string value;
      Parameter *param;
      param = new ScalarParam<std::string>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamString", in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testCArrayParamIntWrite() 
   {
      printMethod(TEST_FUNC);
      int value[3];
      value[0] = 3;
      value[1] = 34;
      value[2] = 8;
      Parameter *param;
      param = new CArrayParam<int>("MyLabel", value, 3);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testCArrayParamIntRead() {
      printMethod(TEST_FUNC);
      int value[3];
      Parameter *param;
      param = new CArrayParam<int>("MyLabel", value, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamInt", in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testCArrayParamDoubleWrite() {
      printMethod(TEST_FUNC);
      double value[3];
      value[0] = 3.0;
      value[1] = 34.7;
      value[2] = 8.97296;
      Parameter *param;
      param = new CArrayParam<double>("MyLabel", value, 3);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testCArrayParamDoubleRead() {
      printMethod(TEST_FUNC);
      double value[3];
      Parameter *param;
      param = new CArrayParam<double>("MyLabel", value, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testFArrayParamIntWrite() 
   {
      printMethod(TEST_FUNC);
      FArray<int, 3> value;
      value[0] = 3;
      value[1] = 34;
      value[2] = 8;
      Parameter *param;
      param = new FArrayParam<int, 3>("MyLabel", value);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testFArrayParamIntRead() {
      printMethod(TEST_FUNC);
      FArray<int, 3> value;
      Parameter *param;
      param = new FArrayParam<int, 3>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ArrayParamInt", in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testFArrayParamDoubleWrite() {
      printMethod(TEST_FUNC);
      FArray<double,3> value;
      value[0] = 3.0;
      value[1] = 34.7;
      value[2] = 8.97296;
      Parameter *param;
      param = new FArrayParam<double, 3>("MyLabel", value);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testFArrayParamDoubleRead() {
      printMethod(TEST_FUNC);
      FArray<double,3> value;
      Parameter *param;
      param = new FArrayParam<double, 3>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testCArray2DParamDoubleWrite() {
      printMethod(TEST_FUNC);
      double value[2][2];
      value[0][0] = 3.0;
      value[0][1] = 34.7;
      value[1][0] = 8.97296;
      value[1][1] = 27.54;
      Parameter *param;
      param = new CArray2DParam<double>("MyLabel", &value[0][0], 2, 2, 2);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testDMatrixParamDoubleWrite() {
      printMethod(TEST_FUNC);
      DMatrix<double> value;
      value.allocate(2, 2);
      value(0, 0) = 3.0;
      value(0, 1) = 34.7;
      value(1, 0) = 8.97296;
      value(1, 1) = 27.54;
      Parameter *param;
      param = new DMatrixParam<double>("MyLabel", value, 2, 2);
      if (verbose() > -3) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testDMatrixParamDoubleRead() {
      printMethod(TEST_FUNC);
      DMatrix<double> value;
      value.allocate(2, 2);
      Parameter *param;
      param = new DMatrixParam<double>("MyLabel", value, 2, 2);
      std::ifstream in;
      openInputFile("in/MatrixParamDouble", in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

};

TEST_BEGIN(ParamTest)
TEST_ADD(ParamTest, testParamIntConstructor)
TEST_ADD(ParamTest, testParamIntWrite)
TEST_ADD(ParamTest, testParamIntRead)
TEST_ADD(ParamTest, testParamDoubleWrite)
TEST_ADD(ParamTest, testParamDoubleRead)
TEST_ADD(ParamTest, testParamStringWrite)
TEST_ADD(ParamTest, testParamStringRead)
TEST_ADD(ParamTest, testCArrayParamIntWrite)
TEST_ADD(ParamTest, testCArrayParamIntRead)
TEST_ADD(ParamTest, testCArrayParamDoubleWrite)
TEST_ADD(ParamTest, testCArrayParamDoubleRead)
TEST_ADD(ParamTest, testFArrayParamIntWrite)
TEST_ADD(ParamTest, testFArrayParamIntRead)
TEST_ADD(ParamTest, testFArrayParamDoubleWrite)
TEST_ADD(ParamTest, testFArrayParamDoubleRead)
TEST_ADD(ParamTest, testCArray2DParamDoubleWrite)
TEST_ADD(ParamTest, testDMatrixParamDoubleWrite)
TEST_ADD(ParamTest, testDMatrixParamDoubleRead)
TEST_END(ParamTest)

#endif
