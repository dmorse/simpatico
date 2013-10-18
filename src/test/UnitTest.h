#ifndef UNIT_TEST_H
#define UNIT_TEST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TestException.h"

#include <fstream>
#include <string>
#ifdef TEST_MPI
#include <mpi.h>
#endif

/**
* UnitTest is a base class for classes that define unit tests.
*
* Each subclass of UnitTest should define one or more test
* methods. Each test method must be a zero parameter function
* that returns void. Test methods may be given arbitrary names.
* Individual test methods should use the preprocessor macro
* TEST_ASSERT(expression) defined in TextException.h to assert 
* the truth of logical expressions.
*
* The test methods defined by a UnitTest are run by an 
* associated subclass of TestRunner. Each test method of a 
* UnitTest must be added to the associated TestRunner.  The 
* run() method of a TestRunner calls all of the associated 
* test methods in the order in which they were added, and 
* counts the number of successful and failed tests. 
*
* The TestRunner associated with a single UnitTest is defined 
* by a class template UnitTestRunner, which takes a UnitTest 
* subclass as a template argument. For example, the TestRunner 
* associated with a UnitTest subclass named TestA is a template 
* instantiation UnitTestRunner<TestA>.
*
* Preprocessor macros defined in the file UnitTestRunner.h 
* should be used to create the boiler-plate code necessary 
* to define a unit test runner and to add test methods to 
* it. 
*/
class UnitTest
{

public:

   /**
   * Constructor
   */
   UnitTest();

   /**
   * Destructor.
   */
   virtual ~UnitTest();

   /**
   * Set up before each test method (empty default implementation).
   */
   virtual void setUp();

   /**
   * Tear down after each test method (empty default implementation).
   */
   virtual void tearDown();

   /**
   * Set verbosity level.
   *
   * \param verbose verbosity level (0 = silent).
   */
   void setVerbose(int verbose);

   /**
   * Set file prefix.
   */
   void setFilePrefix(const std::string& prefix);

   /**
   * Get file prefix string
   */
   const std::string& filePrefix();

   /**
   * Should this processor read and write to file?
   */
   bool isIoProcessor() const;

   #ifdef TEST_MPI
   /**
   * Set the MPI communicator.
   *
   * \param communicator MPI communicator object
   */
   void setCommunicator(MPI::Intracomm& communicator);

   /**
   * Return rank of this processor in communicator.
   */
   int mpiRank();

   /**
   * Does this test have a communicator?
   */
   bool hasCommunicator();

   /**
   * Return the communicator by reference.
   */
   MPI::Intracomm& communicator();
   #endif

protected:

   /**
   * Write name of a class method, iff ioProcessor.
   *
   * \param methodName name of class test method
   */
   void printMethod(const char* methodName);

   /**
   * Write carriage return, iff isIoProcessor.
   */
   void printEndl();

   /**
   * Print a line of hashes, iff isIoProcessor.
   */
   virtual void endMarker();

   /**
   * Open input file.
   *
   * This function adds the filePrefix before the name parameter.
   * It does not check if this node isIoProcessor.
   *
   * \param name base file name (added to filePrefix).
   * \param in input file (opened on return).
   */
   void openInputFile(const std::string& name, std::ifstream& in) const;

   /**
   * Open output file stream.
   *
   * This function adds the filePrefix before the name parameter.
   * It does not check if this node isIoProcessor.
   *
   * \param name base file name (added to filePrefix).
   * \param out  output file (opened on return).
   */
   void openOutputFile(const std::string& name, std::ofstream& out) const;

   /**
   * Return integer verbosity level  (0 == silent).
   */
   int verbose() const;
 
   /**
   * Return true if two integers are equal.
   */
   static bool eq(int s1, int s2);

   /**
   * Return true if two double precision floats are equal.
   */
   static bool eq(double s1, double s2);

private:

   /// Prefix string prepended to file names
   std::string  filePrefix_;

   #ifdef TEST_MPI
   /// Communicator for MPI jobs.
   MPI::Intracomm* communicatorPtr_;

   /// Mpi rank of this processor within communicator.
   int  mpiRank_;
   #endif

   /// Verbosity index
   int  verbose_;

   /// Is this an IoProcessor?
   bool isIoProcessor_;

};

/*
* Constructor
*/
UnitTest::UnitTest() :
   #ifdef TEST_MPI      
   communicatorPtr_(0),
   mpiRank_(-1),
   #endif
   verbose_(0),
   isIoProcessor_(true)
{
   #ifdef TEST_MPI
   // Set the communicator to COMM_WORLD by default.
   setCommunicator(MPI::COMM_WORLD);
   #endif
}

/*
* Destructor.
*/
UnitTest::~UnitTest() 
{}

/*
* Set up before each test method (empty default implementation).
*/
void UnitTest::setUp()
{}

/*
* Tear down after each test method (empty default implementation).
*/
void UnitTest::tearDown()
{}

/*
* Set verbosity level.
*
* \param verbose verbosity level (0 = silent).
*/
void UnitTest::setVerbose(int verbose)
{  verbose_ = verbose; }

/*
* Set file prefix.
*/
void UnitTest::setFilePrefix(const std::string& prefix)
{  filePrefix_  = prefix; }

/*
* Get file prefix string
*/
const std::string& UnitTest::filePrefix()
{  return filePrefix_; }

/*
* Should this processor read and write to file?
*/
bool UnitTest::isIoProcessor() const
{  return isIoProcessor_; } 

#ifdef TEST_MPI
/*
* Set the MPI communicator.
*/
void UnitTest::setCommunicator(MPI::Intracomm& communicator)
{  
   communicatorPtr_ = &communicator; 
   mpiRank_ = communicator.Get_rank();
   if (mpiRank_ == 0) {
      isIoProcessor_ = true; 
   } else {
      isIoProcessor_ = false; 
   }
}

/*
* Return rank of this processor in communicator.
*/
int UnitTest::mpiRank()
{  return mpiRank_; } 

/*
* Does this test have a communicator?
*/
bool UnitTest::hasCommunicator()
{  return bool(communicatorPtr_ != 0); }

/*
* Return the communicator by reference.
*/
MPI::Intracomm& UnitTest::communicator()
{  return *communicatorPtr_; }
#endif

/*
* Write name of a class method, iff ioProcessor.
*/
void UnitTest::printMethod(const char* methodName)
{  if (isIoProcessor()) {
      std::cout << std::endl;
      std::cout << std::string(methodName); 
   }
}

/*
* Write carriage return, iff isIoProcessor.
*/
void UnitTest::printEndl()
{  if (isIoProcessor()) std::cout << std::endl; } 

/*
* Print a line of hashes, iff isIoProcessor.
*/
void UnitTest::endMarker()
{
   if (isIoProcessor()) {
      std::cout << std::endl;
      std::cout << "----------------------------------------------------";
      std::cout << std::endl << std::endl;
   }
}

/*
* Open input file.
*/
void 
UnitTest::openInputFile(const std::string& name, std::ifstream& in) 
const
{   
   std::string filename = filePrefix_;
   filename += name;
   in.open(filename.c_str());
   if (in.fail()) {
      std::cout << std::endl;
      std::cout << "Failure to open input file " 
                << filename << std::endl;
      TEST_THROW("Failure to open file");
   }
}

/*
* Open output file stream.
*/
void 
UnitTest::openOutputFile(const std::string& name, std::ofstream& out) 
const
{   
   std::string filename = filePrefix_;
   filename += name;
   out.open(filename.c_str());
   if (out.fail()) {
      std::cout << std::endl;
      std::cout << "Failure to open output file " 
                << filename << std::endl;
      TEST_THROW("Failure to open file");
   }
}

/*
* Return integer verbosity level  (0 == silent).
*/
inline int UnitTest::verbose() const
{  return verbose_; }
 
/*
* Return true if two integers are equal (static).
*/
inline bool UnitTest::eq(int s1, int s2)
{  return (s1 == s2); }

/*
* Return true if two double precision floats are equal (static).
*/
bool UnitTest::eq(double s1, double s2)
{
   double epsilon = 1.0E-10; 
   return ( fabs(s1-s2) < epsilon ); 
}

#endif
