#ifndef TEST_EXCEPTION_H
#define TEST_EXCEPTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

/**
* An exception thrown by a failed unit test.
*/
class TestException 
{

public:

   /**
   * Default constructor.
   */
   TestException() 
    : message_()
   {}

   /**
   * Constructor for throwing.
   *
   * Constructs error message that includes file and line number. Values 
   * of the file and line parameters should be given by the built-in macros
   * __FILE__ and __LINE__, respectively, in the calling function. A typical
   * call of the constructor is thus of the form:
   * \code
   *    throw TestException("MyClass::myFunction", "A terrible thing happened!",
   *                    __FILE__, __LINE__ );
   * \endcode
   *
   * \param function name of the function from which the TestException was thrown
   * \param message  message describing the nature of the error
   * \param file     name of the file from which the TestException was thrown
   * \param line     line number in file
   */
   TestException(const char *function, const char *message, 
                 const char *file, int line) 
   {
      message_ = "   Function: ";
      message_.append(function);
      message_.append("\n");
      message_.append("   Message : ");
      message_.append(message);
      message_.append("\n");
      message_.append("   File    : ");
      message_.append(file);
      message_.append("\n");
      message_.append("   Line    : ");
      std::ostringstream s;
      s << line;
      message_.append(s.str());
      message_.append("\n");
   }
   
   /**
   * Constructor without function name parameter.
   *
   * \param message  message describing the nature of the error
   * \param file     name of the file from which the TestException was thrown
   * \param line     line number in file
   */
   TestException(const char *message, const char *file, int line) 
   {
      message_ = "   Message : ";
      message_.append(message);
      message_.append("\n");
      message_.append("   File    : ");
      message_.append(file);
      message_.append("\n");
      message_.append("   Line    : ");
      std::ostringstream s;
      s << line;
      message_.append(s.str());
      message_.append("\n");
   }
   
   /**
   * Destructor.
   */
   ~TestException()
   {}
   
   /**
   * Write error message to output stream.
   *
   * \param out output stream
   */
   void write(std::ostream &out)
   {  out << message_; }
   
   /**
   * Return the error message.
   */
   const std::string& message()
   {  return message_; }


protected:

   std::string message_;  ///< Error message

};

/**
* Macro for the name of the current function (compiler dependent).
*/
#define TEST_FUNC __PRETTY_FUNCTION__

/**
* Assert macro, throws a TestException if assertion fails.
*/
#ifdef  TEST_FUNC
  #define TEST_ASSERT(expr) \
  if (!(expr)) throw TestException(TEST_FUNC, #expr , __FILE__, __LINE__)
#else
  #define TEST_ASSERT(expr) \
  if (!(expr)) throw TestException( #expr , __FILE__, __LINE__)
#endif

/**
* Throw a TestException, with a message.
*/
#ifdef  TEST_FUNC
  #define TEST_THROW(msg) \
  throw TestException(TEST_FUNC, #msg , __FILE__, __LINE__)
#else
  #define TEST_THROW(msg) \
  throw TestException( #msg , __FILE__, __LINE__)
#endif

#endif
