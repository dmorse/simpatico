#ifndef UTIL_PARAMETER_H
#define UTIL_PARAMETER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>
#include <util/param/Label.h>

class ParamTest;

namespace Util
{

   /**
   * A single variable in a parameter file.
   *
   * Parameter is a base class for objects that read and write the value of
   * a single C++ variable from or to file. Each Parameter has a Label that
   * contains a std::string label that identifies a variable in a parameter
   * file. Instances of subclass also contain a pointer to a C++ variable
   * of the appropriate type.
   *
   * \ingroup Param_Module
   */
   class Parameter : public ParamComponent
   {

   public:

      /// Width of data field for a scalar variable (numerical or string).
      static const int Width = 20;

      /// Precision for io of floating point data field.
      static const int Precision = 12;

      /**
      * Constructor.
      *
      * \param label label string that precedes value in file format.
      */
      Parameter(const char *label);

      /// Destructor.
      virtual ~Parameter();

      /// Return label string.
      std::string label();

   protected:

      /// Label object that contains parameter label string.
      Label label_;

   //friends:

      friend class ::ParamTest;

   };

}
#endif
