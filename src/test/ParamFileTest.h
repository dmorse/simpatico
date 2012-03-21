#ifndef PARAM_FILE_TEST_H
#define PARAM_FILE_TEST_H

#include "UnitTest.h"

#include <string>
#include <iostream>
#include <fstream>

template <class T>
class ParamFileTest : public UnitTest 
{

public:

   /**
   * Constructor.
   */
   ParamFileTest()
    : UnitTest()
   {}

   /**
   * Close the input file.
   */
   virtual void tearDown()
   {  closeFile(); }

   /**
   * Open the input file.
   */
   void openFile(const char *fileName)
   { 
      if (isIoProcessor()) {
         openInputFile(std::string(fileName), file_);
      }
   }

   /**
   * Close the input file.
   */
   void closeFile()
   {
      if (file_.is_open())
         file_.close(); 
   }

   /**
   * Returns associated object by reference.
   */
   T& object()
   {  return object_; }

   /**
   * Returns input file by reference.
   */
   std::ifstream& file()
   {  return file_; }

private:

   T              object_;
   std::ifstream  file_;

};

#endif
