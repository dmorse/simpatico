#ifndef UNIT_TEST_RUNNER_H
#define UNIT_TEST_RUNNER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TestRunner.h"
#include "TestException.h"
#include <vector>

/**
* Template for a TestRunner that runs an associated UnitTest.
*
* A instance of UnitTestRunner<MyTest> holds an array of pointers to
* the test methods of a class MyTest that is derived from UnitTest.
* Each of the test methods must return void and take no parameters. The
* addTestMethod() method adds a pointer to such a method to this array. 
* The run() method runs all of the registered test methods in sequence.
*
* To run a set of unit tests one must: 
*
*  - Define a subclass of UnitTest,
*  - Define an associated subclass of UnitTestRunner, 
*  - Construct a UnitTestRunner object and call run().
*
* The code to define the UnitTestRunner class may be simplified with
* a set of preprocessor macros. 
*
* Here is an example of the code to run all the methods of a UnitTest,
* with no preprocessor macros:
* \code
*
* // Define a UnitTest class
* class MyTest : public UnitTest {
* public:
*
*    test1()
*    { .... }
*
*    test2
*    { ...  }
*
* };
*
* // Define a UnitTestRunner associated with MyTest
* class MyTest_Runner : public UnitTestRunner<MyTest> {
* public:
*
*    MyTest_Runner(){
*       addTestMethod(&MyTest::test1);
*       addTestMethod(&MyTest::test2);
*    }
*
* }
*
* // Run the tests.
* MyTest_Runner runner;
* runner.run();
*
* \endcode
* The following series of preprocessor macros may be used to generate 
* the definition of the MyTest_Runner class in the above example, and 
* to create an instance of this class:
* \code
* 
* TEST_BEGIN(MyTest)
* TEST_ADD(MyTest, test1)
* TEST_ADD(MyTest, test2)
* TEST_END(MyTest)
*
* TEST_RUNNER(MyTest) runner;
* runner.run();
*
* \endcode
* The macro TEST_BEGIN(TestClass) generates the beginning of the above
* class definition for a new subclass of UnitTestRunner<TestClass>. The
* TEST_ADD(TestClass, Method) adds a specified method of the associated
* class TestClass to the constructor of the new UnitTestRunner class. 
* The TEST_END macro closes the class definition. The name of the
* resulting UnitTestRunner class is given by the preprocessor macro 
* TEST_RUNNER(TestClass). This macro may thus be used as a class name 
* to instantiate objects of this class. By convention, the macro
* TEST_RUNNER(TestClass) expands to the name TestClass_Runner, in which
* a standard suffix "_Runner" is added to the name of the UnitTest 
* class.
*/
template <class UnitTestClass>
class UnitTestRunner : public TestRunner
{

public:

   using TestRunner::nFailure;
   using TestRunner::isIoProcessor;

   /**
   * Pointer to a test method of the associated UnitTest class.
   */
   typedef void (UnitTestClass::*MethodPtr)();

   // Default constructor.
   UnitTestRunner()
    : TestRunner()
   {
      #ifdef TEST_MPI
      if (isIoProcessor()) {
         results_.reserve(mpiSize());
         for (int i=0; i < mpiSize(); ++i) {
            results_.push_back(false);
         }
      }
      #endif
   }

   // Default destructor.

   /**
   * Register a test method of associated test case.
   */
   void addTestMethod(MethodPtr methodPtr)
   {  methodPtrs_.push_back(methodPtr); }

   /**
   * Return number of registered test methods.
   */
   int nTestMethod()
   {  return methodPtrs_.size(); } 

   /**
   * Run test method number i.
   *
   * \param i index of test method
   */
   void method(unsigned int i)
   {
      UnitTestClass testCase;
      #ifdef TEST_MPI
      TestException exception;
      int result;
      #endif

      testCase.setFilePrefix(filePrefix());
      try {
         testCase.setUp();
         (testCase.*methodPtrs_[i])();
         #ifndef TEST_MPI
         if (testCase.isIoProcessor()) {
            recordSuccess();
         }
         #else
         result = 1;
         #endif
      } catch (TestException &e) {
         #ifndef TEST_MPI
         std::cout << std::endl;
         std::cout << " Failure " << std::endl << std::endl;
         std::cout << e.message() << std::endl;
         recordFailure();
         #else
         result = 0;
         exception = e;
         #endif
      }
      testCase.tearDown();

      #ifndef TEST_MPI
      std::cout << ".";
      #else
      MPI::COMM_WORLD.Barrier();
      if (mpiRank() == 0) {
         results_[0] = result;
         if (results_[0] == 0) {
            std::cout << std::endl;
            std::cout << " Failure  on Processor 0" 
                      << std::endl << std::endl;
            std::cout << exception.message() << std::endl;
            std::cout.flush();
         }
         for (int i=1; i < mpiSize(); ++i) {
            MPI::COMM_WORLD.Recv(&(results_[i]), 1, MPI_INT, i, i);
            if (results_[i] == 0) {
               result = 0;
               MPI::COMM_WORLD.Send(&(results_[i]), 1, MPI_INT, i, mpiSize() + i);
               MPI::COMM_WORLD.Recv(&(results_[i]), 1, MPI_INT, i, 2*mpiSize() + i);
            }
         }
         if (result) {
            recordSuccess();
         } else {
            recordFailure();
         }
         std::cout << ".";
         std::cout.flush();
      } else {
         MPI::COMM_WORLD.Send(&result, 1, MPI_INT, 0, mpiRank());
         if (result == 0) {
            MPI::COMM_WORLD.Recv(&result, 1, MPI_INT, 0, 
                                 mpiSize() + mpiRank());
            std::cout.flush();
            std::cout << std::endl;
            std::cout << " Failure  on Processor " << mpiRank() 
                      << std::endl << std::endl;
            std::cout << exception.message() << std::endl;
            std::cout.flush();
            MPI::COMM_WORLD.Send(&result, 1, MPI_INT, 0, 
                                 2*mpiSize() + mpiRank());
         }
      }
      MPI::COMM_WORLD.Barrier();
      #endif
   }

   /**
   * Run all registered test methods sequentially.
   */
   virtual int run()
   {
      for (unsigned int i = 0; i < methodPtrs_.size(); ++i) {
         method(i); 
      }
      report();
      return nFailure();
   }

private:

   std::vector<MethodPtr> methodPtrs_;

   #ifdef TEST_MPI
   std::vector<int> results_;
   #endif 

};

/**
* Macro for name of the UnitTestRunner class associated with UnitTestClass.
*/
#define TEST_RUNNER(UnitTestClass) UnitTestClass##_Runner

/**
* Begin definition of class TEST_RUNNER(UnitTestClass).
*
* This macro generates code to open the class definition, and to open
* the definition of a default constructor.
*/
#define TEST_BEGIN(UnitTestClass) \
   class TEST_RUNNER(UnitTestClass) \
    : public UnitTestRunner<UnitTestClass> \
   { public: TEST_RUNNER(UnitTestClass)() { 

/**
* Macro to add a test method to TEST_RUNNER(UnitTestClass).
*/
#define TEST_ADD(UnitTestClass, Method) \
   addTestMethod(&UnitTestClass::Method);

/**
* Macro to end definition of a class TEST_RUNNER(UnitTestClass).
*
* This macro ends both the constructor and class definition.
*/
#define TEST_END(UnitTestClass) } }; 

#endif
