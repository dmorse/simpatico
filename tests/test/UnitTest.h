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
* that returns void, and may be given any informative name.
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
   UnitTest() :
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

   /**
   * Destructor.
   */
   virtual ~UnitTest() 
   {}

   /**
   * Set up before each test method (empty default implementation).
   */
   virtual void setUp()
   {}

   /**
   * Tear down after each test method (empty default implementation).
   */
   virtual void tearDown()
   {}

   /**
   * Set verbosity level.
   *
   * \param verbose verbosity level (0 = silent).
   */
   void setVerbose(int verbose)
   {  verbose_ = verbose; }

   /**
   * Set file prefix.
   */
   void setFilePrefix(const std::string& prefix)
   {  filePrefix_  = prefix; }

   /**
   * Get file prefix string
   */
   const std::string& filePrefix()
   {  return filePrefix_; }

   /**
   * Should this processor read and write to file?
   */
   bool isIoProcessor() const
   {  return isIoProcessor_; } 

   #ifdef TEST_MPI
   void setCommunicator(MPI::Intracomm& communicator)
   {  
      communicatorPtr_ = &communicator; 
      mpiRank_ = communicator.Get_rank();
      if (mpiRank_ == 0) {
         isIoProcessor_ = true; 
      } else {
         isIoProcessor_ = false; 
      }
   }

   /**
   * Return rank of this processor in communicator.
   */
   int mpiRank()
   {  return mpiRank_; } 

   /**
   * Does this test have a communicator?
   */
   bool hasCommunicator()
   {  return bool(communicatorPtr_ != 0); }

   /**
   * Return the communicator by reference.
   */
   MPI::Intracomm& communicator()
   {  return *communicatorPtr_; }
   #endif

protected:

   /**
   * Write name of a class method, iff ioProcessor.
   */
   void printMethod(const char* methodName)
   {  if (isIoProcessor()) {
         std::cout << std::endl;
         std::cout << std::string(methodName); 
      }
   }

   /**
   * Write carriage return, iff isIoProcessor.
   */
   void printEndl()
   {  if (isIoProcessor()) std::cout << std::endl; } 

   /**
   * Print a line of hashes, iff isIoProcessor.
   */
   virtual void endMarker()
   {
      if (isIoProcessor()) {
         std::cout << std::endl;
         std::cout << "----------------------------------------------------";
         std::cout << std::endl << std::endl;
      }
   }

   /**
   * Open input file.
   *
   * This function adds the filePrefix before the name parameter.
   * It does not check if this node isIoProcessor.
   *
   * \param name base file name (added to filePrefix).
   * \param in input file (opened on return).
   */
   void openInputFile(const std::string& name, std::ifstream& in) const
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

   /**
   * Open output file stream.
   *
   * This function adds the filePrefix before the name parameter.
   * It does not check if this node isIoProcessor.
   *
   * \param name base file name (added to filePrefix).
   * \param out  output file (opened on return).
   */
   void openOutputFile(const std::string& name, std::ofstream& out) const
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

   /**
   * Return integer verbosity level  (0 == silent).
   */
   int verbose() const
   {  return verbose_; }
 
   /**
   * Return true if two integers are equal.
   */
   static bool eq(int s1, int s2)
   {  return (s1 == s2); }

   /**
   * Return true if two double precision floats are equal.
   */
   static bool eq(double s1, double s2)
   {
      double epsilon = 1.0E-10; 
      return ( fabs(s1-s2) < epsilon ); 
   }

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
#endif
