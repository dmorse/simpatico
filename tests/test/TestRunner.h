#ifndef TEST_RUNNER_H
#define TEST_RUNNER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>

#ifdef TEST_MPI
#include <mpi.h>
#endif

/**
* Abstract base class for classes that run tests.
*
* TestRunner is an abstract base class with two types of subclass:
* The UnitTestRunner class template defines a TestRunner that runs
* the tests for an associated UnitTest. A CompositeTestRunner runs 
* the tests for a sequence of other TestRunner objects, each of 
* which can be a UnitTestRunner or another CompositeTestRunner.
*
* The pure virtual run() method of a TestRunner must run all of 
* the associated test methods, and records the number nSuccess() 
* of tests that succeed and the number nFailure() that fail. A 
* test fails if it throws a TestException. Must test methods use
* the TEST_ASSERT(expr) macro to throw a TestException if the
* logical expression expr is false.
*/
class TestRunner
{

public:

   /**
   * Constructor.
   */
   TestRunner()
    : parentPtr_(0),
      nSuccess_(0),
      nFailure_(0),
      isIoProcessor_(0)
      #ifdef TEST_MPI
      ,mpiRank_(0),
      mpiSize_(0)
      #endif
   {
       #ifndef TEST_MPI
          isIoProcessor_ = true; 
       #else
          mpiRank_ = MPI::COMM_WORLD.Get_rank();
          mpiSize_ = MPI::COMM_WORLD.Get_size();
          if (mpiRank_ == 0) {
             isIoProcessor_ = true; 
          } else {
             isIoProcessor_ = false; 
          }
       #endif
   }

   /**
   * Destructor.
   */
   virtual ~TestRunner()
   {}

   /**
   * Run all tests.
   *
   * \return number of failures.
   */
   virtual int run() = 0;

   /**
   * Increment counter for failed tests, and that of parent (if any).
   */
   void recordFailure() 
   {
      if (isIoProcessor()) {
         ++nFailure_;
         if (hasParent()) {
            parent().recordFailure();
         }
      }
   }

   /**
   * Increment counter for successful tests, and that of parent (if any).
   */
   void recordSuccess()
   {
      if (isIoProcessor()) {
         ++nSuccess_;
         if (hasParent()) {
            parent().recordSuccess();
         }
      }
   }

   /**
   * Set another TestRunner as the parent.
   */
   void setParent(TestRunner& parent)
   {  parentPtr_ = &parent; }

   /**
   * Return the parent object, if any.
   */
   TestRunner& parent()
   {  return *parentPtr_; }

   /**
   * Does this object have a parent?
   */
   bool hasParent() const
   {  return (parentPtr_ != 0); }

   /**
   * Return number of successful tests run.
   */
   int nSuccess() const
   {  return nSuccess_; }

   /**
   * Return number of failed tests run.
   */
   int nFailure() const
   {  return nFailure_; }

   /**
   * If this object has no parent, report success and failure counters.
   */
   void report() const
   { 
      if (!hasParent() && isIoProcessor()) {
         std::cout << std::endl;
         std::cout << nSuccess_ << "  successful tests  " << std::endl;
         std::cout << nFailure_ << "  failed tests  "    << std::endl;
         std::cout << std::endl;
      }
   }

   bool isIoProcessor() const
   {  return isIoProcessor_; }

   #ifdef TEST_MPI
   int mpiRank() const
   {  return mpiRank_; }

   int mpiSize() const
   {  return mpiSize_; } 
   #endif

   /**
   * Prepend argument prefix to existing filePrefix.
   */
   virtual void addFilePrefix(const std::string& prefix) 
   {
      std::string newPrefix = prefix;
      newPrefix += filePrefix_;
      filePrefix_ = newPrefix;
   }

   /**
   * Return file prefix by const reference.
   */
   const std::string& filePrefix() const
   {  return filePrefix_; }

protected:

   /// Prefix added to file names
   std::string  filePrefix_;

private:

   /// Pointer to a parent TestRunner (if any).
   TestRunner* parentPtr_;

   /// Total number of successful tests run.
   int  nSuccess_;

   /// Total number of failed tests run.
   int  nFailure_;

   /// Can this processor input and output data?
   /// This is always true when TEST_MPI is not defined.
   bool  isIoProcessor_;

   #ifdef TEST_MPI
   /// Rank of this processor within an MPI job.
   int  mpiRank_;

   /// Size of associated MPI communicator.
   int  mpiSize_;
   #endif

};
#endif
