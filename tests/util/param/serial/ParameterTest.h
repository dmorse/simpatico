#ifndef PARAMETER_TEST_H
#define PARAMETER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/param/Parameter.h>
#include <util/param/ScalarParam.h>
#include <util/param/CArrayParam.h>
#include <util/param/DArrayParam.h>
#include <util/param/FArrayParam.h>
#include <util/param/CArray2DParam.h>
#include <util/param/DMatrixParam.h>

using namespace Util;

class ParameterTest : public UnitTest 
{

public:

   void setUp()
   { 
      setVerbose(2); 
      ParamComponent::setEcho(true);
   }

   void tearDown()
   {
      ParamComponent::setEcho(false);
   }

   void testParamIntConstructor1() {
      printMethod(TEST_FUNC);
      int        value = 4;
      Parameter* param;
      param = new ScalarParam<int>("MyLabel", value);
      delete param;
   }

   void testParamIntConstructor2() {
      printMethod(TEST_FUNC);
      int        value = 4;
      Parameter* param;
      param = new ScalarParam<int>("MyLabel", value, false);
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

   void testParamIntRead1() 
   {
      printMethod(TEST_FUNC);
      int        value;
      Parameter *param;
      param  = new ScalarParam<int>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      in.close();
      TEST_ASSERT(Label::isClear());
      TEST_ASSERT(param->label() == "MyLabel");
      TEST_ASSERT(value == 36);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testParamIntRead2() 
   {
      /// Test absent optional parameter
      printMethod(TEST_FUNC);
      int empty;
      int value;
      Parameter* absent = new ScalarParam<int>("Absent", empty, false);
      Parameter* param  = new ScalarParam<int>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absent->readParam(in);
      TEST_ASSERT(!Label::isClear());
      param->readParam(in);
      TEST_ASSERT(Label::isClear());
      TEST_ASSERT(param->label() == "MyLabel");
      //TEST_ASSERT(value == 36);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      in.close();
      delete absent;
      delete param;
   }

   void testParamIntReadSaveLoad1() 
   {
      printMethod(TEST_FUNC);
      int        value;
      Parameter *param;
      param  = new ScalarParam<int>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      in.close();
      TEST_ASSERT(Label::isClear());
      TEST_ASSERT(param->label() == "MyLabel");
      TEST_ASSERT(value == 36);

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      param->save(oar);
      oar.file().close();

      // Load from archive
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      int value2;
      Parameter* param2  = new ScalarParam<int>("MyLabel", value2);
      param2->load(iar);
      iar.file().close();
      TEST_ASSERT(param2->label() == "MyLabel");
      TEST_ASSERT(value2 == 36);

      if (verbose() > 0) {
         printEndl();
         param2->writeParam(std::cout);
      }

      delete param;
      delete param2;
   }

   void testParamIntReadSaveLoad2() 
   {
      printMethod(TEST_FUNC);
      int empty;
      int value;
      Parameter* absent = new ScalarParam<int>("Wrong", empty, false);
      Parameter* param  = new ScalarParam<int>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absent->readParam(in);
      TEST_ASSERT(!Label::isClear());
      param->readParam(in);
      in.close();
      TEST_ASSERT(Label::isClear());
      TEST_ASSERT(param->label() == "MyLabel");
      TEST_ASSERT(value == 36);

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      absent->save(oar);
      param->save(oar);
      oar.file().close();

      // Load from archive
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      if (ParamComponent::echo()) std::cout << std::endl;
      int empty2;
      int value2;
      Parameter* absent2 = new ScalarParam<int>("Wrong", empty2, false);
      Parameter* param2 = new ScalarParam<int>("MyLabel", value2);
      absent2->load(iar);
      param2->load(iar);
      iar.file().close();
      TEST_ASSERT(param2->label() == "MyLabel");
      TEST_ASSERT(value2 == 36);

      if (verbose() > 0) {
         printEndl();
         param2->writeParam(std::cout);
      }

      delete param;
      delete absent;
      delete param2;
      delete absent2;
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
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      in.close();
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

   void testParamStringRead1() {
      printMethod(TEST_FUNC);
      std::string value;
      Parameter *param;
      param = new ScalarParam<std::string>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamString", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testParamStringRead2() {
      printMethod(TEST_FUNC);
      int empty;
      std::string value;
      Parameter* absent = new ScalarParam<int>("Wrong", empty, false);
      Parameter* param = new ScalarParam<std::string>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamString", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absent->readParam(in);
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete absent;
      delete param;
   }

   void testParamStringReadSaveLoad2() 
   {
      printMethod(TEST_FUNC);
      int empty;
      std::string value;
      Parameter* absent = new ScalarParam<int>("Wrong", empty, false);
      Parameter* param  = new ScalarParam<std::string>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ScalarParamString", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absent->readParam(in);
      TEST_ASSERT(!Label::isClear());
      param->readParam(in);
      in.close();
      TEST_ASSERT(Label::isClear());
      TEST_ASSERT(param->label() == "MyLabel");
      //TEST_ASSERT(value == 36);

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      absent->save(oar);
      param->save(oar);
      oar.file().close();

      // Load from archive
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      int empty2;
      std::string value2;
      Parameter* absent2 = new ScalarParam<int>("Wrong", empty2, false);
      Parameter* param2 = new ScalarParam<std::string>("MyLabel", value2);
      absent2->load(iar);
      param2->load(iar);
      iar.file().close();
      TEST_ASSERT(param2->label() == "MyLabel");
      //TEST_ASSERT(value2 == 36);

      if (verbose() > 0) {
         printEndl();
         param2->writeParam(std::cout);
      }

      delete param;
      delete absent;
      delete param2;
      delete absent2;
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
      if (ParamComponent::echo()) std::cout << std::endl;
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
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testDArrayParamIntWrite() 
   {
      printMethod(TEST_FUNC);
      DArray<int> value;
      value.allocate(3);
      value[0] = 3;
      value[1] = 34;
      value[2] = 8;
      Parameter* param = new DArrayParam<int>("MyLabel", value, 3);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testDArrayParamIntRead() {
      printMethod(TEST_FUNC);
      DArray<int> value;
      value.allocate(3);
      Parameter* param = new DArrayParam<int>("MyLabel", value, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testDArrayParamDoubleWrite() {
      printMethod(TEST_FUNC);
      DArray<double> value;
      value.allocate(3);
      value[0] = 3.0;
      value[1] = 34.7;
      value[2] = 8.97296;
      Parameter *param;
      param = new DArrayParam<double>("MyLabel", value, 3);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testDArrayParamDoubleRead() {
      printMethod(TEST_FUNC);
      DArray<double> value;
      value.allocate(3);
      Parameter *param;
      param = new DArrayParam<double>("MyLabel", value, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testDArrayParamDoubleReadSaveLoad1() {
      printMethod(TEST_FUNC);
      DArray<double> value;
      value.allocate(3);
      Parameter* param  = new DArrayParam<double>("MyLabel", value, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      param->save(oar);
      oar.file().close();

      // Load from archive
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      DArray<double> value2;
      value2.allocate(3);
      Parameter* param2 = new DArrayParam<double>("MyLabel", value2, 3);
      param2->load(iar);
      iar.file().close();
      TEST_ASSERT(param2->label() == "MyLabel");
      //TEST_ASSERT(value2 == 36);

      if (verbose() > 0) {
         printEndl();
         param2->writeParam(std::cout);
      }

      delete param;
      //delete param2;
   }

   void testDArrayParamDoubleReadSaveLoad2() {
      printMethod(TEST_FUNC);
      DArray<double> empty;
      DArray<double> value;
      empty.allocate(3);
      value.allocate(3);
      Parameter* absent = new DArrayParam<double>("Wrong", empty, 3, false);
      Parameter* param  = new DArrayParam<double>("MyLabel", value, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absent->readParam(in);
      param->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      absent->save(oar);
      param->save(oar);
      oar.file().close();

      // Load from archive
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      DArray<double> empty2;
      DArray<double> value2;
      empty2.allocate(3);
      value2.allocate(3);
      Parameter* absent2 = new DArrayParam<double>("Wrong", empty2, 3, false);
      Parameter* param2 = new DArrayParam<double>("MyLabel", value2, 3);
      absent2->load(iar);
      param2->load(iar);
      iar.file().close();

      TEST_ASSERT(param2->label() == "MyLabel");
      //TEST_ASSERT(value2 == 36);

      if (verbose() > 0) {
         printEndl();
         param2->writeParam(std::cout);
      }

      delete absent;
      delete param;
      delete absent2;
      delete param2;
   }

   void testDArrayParamDoubleReadSaveLoad3() {
      printMethod(TEST_FUNC);
      DArray<double> empty;
      DArray<double> value;
      empty.allocate(3);
      value.allocate(3);
      Parameter* absent = new DArrayParam<double>("Wrong", empty, 3, false);
      Parameter* param  = new DArrayParam<double>("MyLabel", value, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absent->readParam(in);
      param->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptional(oar, empty, false);
      Parameter::saveOptional(oar, value, true);
      //absent->save(oar);
      //param->save(oar);
      oar.file().close();

      #if 0
      // Load from archive
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      DArray<double> empty2;
      DArray<double> value2;
      empty2.allocate(3);
      value2.allocate(3);
      Parameter* absent2 = new DArrayParam<double>("Wrong", empty2, 3, false);
      Parameter* param2 = new DArrayParam<double>("MyLabel", value2, 3);
      absent2->load(iar);
      param2->load(iar);
      iar.file().close();

      TEST_ASSERT(param2->label() == "MyLabel");
      //TEST_ASSERT(value2 == 36);

      if (verbose() > 0) {
         printEndl();
         param2->writeParam(std::cout);
      }
      #endif

      delete absent;
      delete param;
      //delete absent2;
      //delete param2;
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
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      in.close();
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
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testFArrayParamDoubleReadSaveLoad1() {
      printMethod(TEST_FUNC);
      FArray<double, 3> value;
      Parameter* param  = new FArrayParam<double, 3>("MyLabel", value);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      param->save(oar);
      oar.file().close();

      // Load from archive
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      FArray<double, 3> value2;
      Parameter* param2 = new FArrayParam<double, 3>("MyLabel", value2);
      param2->load(iar);
      iar.file().close();
      TEST_ASSERT(param2->label() == "MyLabel");
      //TEST_ASSERT(value2 == 36);

      if (verbose() > 0) {
         printEndl();
         param2->writeParam(std::cout);
      }

      delete param;
      //delete param2;
   }

   void testFArrayParamDoubleReadSaveLoad2() {
      printMethod(TEST_FUNC);
      FArray<double, 3> empty;
      FArray<double, 3> value;
      Parameter* absent = new FArrayParam<double, 3>("Wrong", empty, false);
      Parameter* param  = new FArrayParam<double, 3>("MyLabel", value, true);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absent->readParam(in);
      param->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      absent->save(oar);
      param->save(oar);
      oar.file().close();

      // Load from archive
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      FArray<double, 3> empty2;
      FArray<double, 3> value2;
      Parameter* absent2 = new FArrayParam<double, 3>("Wrong", empty2, false);
      Parameter* param2 = new FArrayParam<double, 3>("MyLabel", value2);
      absent2->load(iar);
      param2->load(iar);
      iar.file().close();

      TEST_ASSERT(param2->label() == "MyLabel");
      //TEST_ASSERT(value2 == 36);

      if (verbose() > 0) {
         printEndl();
         param2->writeParam(std::cout);
      }

      delete absent;
      delete param;
      delete absent2;
      delete param2;
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
      if (ParamComponent::echo()) std::cout << std::endl;
      param->readParam(in);
      if (verbose() > 0) {
         printEndl();
         param->writeParam(std::cout);
      }
      delete param;
   }

};

TEST_BEGIN(ParameterTest)
TEST_ADD(ParameterTest, testParamIntConstructor1)
TEST_ADD(ParameterTest, testParamIntConstructor2)
TEST_ADD(ParameterTest, testParamIntWrite)
TEST_ADD(ParameterTest, testParamIntRead1)
TEST_ADD(ParameterTest, testParamIntRead2)
TEST_ADD(ParameterTest, testParamIntReadSaveLoad1)
TEST_ADD(ParameterTest, testParamIntReadSaveLoad2)
TEST_ADD(ParameterTest, testParamDoubleWrite)
TEST_ADD(ParameterTest, testParamDoubleRead)
TEST_ADD(ParameterTest, testParamStringWrite)
TEST_ADD(ParameterTest, testParamStringRead1)
TEST_ADD(ParameterTest, testParamStringRead2)
TEST_ADD(ParameterTest, testParamStringReadSaveLoad2)
TEST_ADD(ParameterTest, testCArrayParamIntWrite)
TEST_ADD(ParameterTest, testCArrayParamIntRead)
TEST_ADD(ParameterTest, testCArrayParamDoubleWrite)
TEST_ADD(ParameterTest, testCArrayParamDoubleRead)
TEST_ADD(ParameterTest, testDArrayParamIntWrite)
TEST_ADD(ParameterTest, testDArrayParamIntRead)
TEST_ADD(ParameterTest, testDArrayParamDoubleWrite)
TEST_ADD(ParameterTest, testDArrayParamDoubleRead)
TEST_ADD(ParameterTest, testDArrayParamDoubleReadSaveLoad1)
TEST_ADD(ParameterTest, testDArrayParamDoubleReadSaveLoad2)
TEST_ADD(ParameterTest, testDArrayParamDoubleReadSaveLoad3)
TEST_ADD(ParameterTest, testFArrayParamIntWrite)
TEST_ADD(ParameterTest, testFArrayParamIntRead)
TEST_ADD(ParameterTest, testFArrayParamDoubleWrite)
TEST_ADD(ParameterTest, testFArrayParamDoubleRead)
TEST_ADD(ParameterTest, testFArrayParamDoubleReadSaveLoad1)
TEST_ADD(ParameterTest, testFArrayParamDoubleReadSaveLoad2)
TEST_ADD(ParameterTest, testCArray2DParamDoubleWrite)
TEST_ADD(ParameterTest, testDMatrixParamDoubleWrite)
TEST_ADD(ParameterTest, testDMatrixParamDoubleRead)
TEST_END(ParameterTest)

#endif
