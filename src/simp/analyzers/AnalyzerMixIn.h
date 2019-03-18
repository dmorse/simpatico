#ifndef SIMP_ANALYZER_MIXIN_H
#define SIMP_ANALYZER_MIXIN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <iostream>
#include <string>

namespace Util { 
   class FileMaster;
}

namespace Simp
{

   using namespace Util;

   /**
   * Base class for Analyzer MixIn classes.
   *
   * An AnalyzerMixIn has an output file and a FileMaster that
   * can be used to open it with the correct prefix.
   *
   * \ingroup Simp_Analyzer_Module
   */
   class AnalyzerMixIn 
   {
   
   public:
   
      /**
      * Constructor.
      *
      * Stores a private pointer to the fileMaster, used by the 
      * openOutputFile function.
      *
      * \param fileMaster associated FileMaster object
      */
      AnalyzerMixIn(FileMaster& fileMaster);
   
      /**
      * Destructor.
      */
      virtual ~AnalyzerMixIn(); 

   protected:

      /**
      * Access output file by reference.
      */
      std::ofstream& outputFile();

      /**
      * Open the output file. 
      *
      * \param filename base file name, without output prefix
      */
      void openOutputFile(std::string filename);

   private:

      /// Output file stream.
      std::ofstream  outputFile_;

      /// Pointer to a FileMaster
      FileMaster* fileMasterPtr_;

   };

   /*
   * Access output file by reference.
   */
   inline
   std::ofstream& AnalyzerMixIn::outputFile()
   {  return outputFile_; }

}
#endif 
