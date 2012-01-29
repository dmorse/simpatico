#ifndef EXCEPTION_H
#define EXCEPTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <string>

namespace Util
{

   /**
   * A user-defined exception.
   *
   * Exceptions are usually thrown using the UTIL_THROW() macro. 
   *
   * \ingroup Util_Module
   */
   class Exception 
   {
   
   public:
   
      /**
      * Constructor. 
      *
      * Constructs error message that includes file and line number. Values 
      * of the file and line parameters should be given by the built-in macros
      * __FILE__ and __LINE__, respectively, in the calling function. A typical
      * call of the constructor is thus of the form:
      * \code
      *    throw Exception("MyClass::myFunction", "A terrible thing happened!",
      *                    __FILE__, __LINE__ );
      * \endcode
      *
      * \param function name of the function from which the Exception was thrown
      * \param message  message describing the nature of the error
      * \param file     name of the file from which the Exception was thrown
      * \param line     line number in file
      * \param echo     if echo, then echo to Log::file() when constructed. 
      */
      Exception(const char *function, const char *message, 
                const char *file, int line, int echo = 1);
   
      /**
      * Constructor without function name parameter.
      *
      * \param message  message describing the nature of the error
      * \param file     name of the file from which the Exception was thrown
      * \param line     line number in file
      * \param echo     if echo, then echo to std out when constructed. 
      */
      Exception(const char *message, const char *file, int line, int echo = 1);

      /**
      * Destructor
      */   
      virtual ~Exception();

      /**
      * Write error message to output stream.
      *
      * \param out output stream
      */
      void write(std::ostream &out);
   
      /**
      * Return the error message.
      */
      const std::string& message();
   
   
   protected:
  
      /// Error message string.
      std::string message_; 
   
   };

}
#endif
