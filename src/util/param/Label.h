#ifndef UTIL_LABEL_H
#define UTIL_LABEL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <string>

namespace Util
{

   /**
   * A label string in a file format.
   *
   * The operator >> for a label checks if the expected label was found.
   *
   * \ingroup Param_Module
   */
   class Label
   {

   public:

      /// Width of label field in file format.
      static const int LabelWidth  = 20;

      /**
      * Constructor.
      *
      * \param label label string that precedes value in file format
      */
      explicit Label(const char* label);

      /**
      * Copy constructor.
      *
      * \param other Label object being copied.
      */
      Label(const Label& other);

      /**
      * Destructor.
      */
      ~Label();

      /**
      * Return label string.
      */
      std::string string() const;

   private:

      std::string label_;

   //friends:

      friend std::istream& operator>>(std::istream& in, Label label);
      friend std::ostream& operator<<(std::ostream& out, Label label);

   };

   /**
   * Extractor for Label.
   *
   * \param in    input stream
   * \param label Label to be read from file
   */ 
   std::istream& operator>>(std::istream& in, Label label);

   /**
   * Inserter for Label.
   *
   * \param out   output stream
   * \param label Label to be written to file
   */ 
   std::ostream& operator<<(std::ostream& out, Label label);

} 
#endif
