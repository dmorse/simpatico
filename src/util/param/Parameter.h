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

namespace Util
{

   /**
   * A single variable in a parameter file.
   *
   * Parameter is a base class for objects that read and write the 
   * value of a single C++ variable from or to file. Each Parameter 
   * has a Label that contains a std::string label that identifies 
   * a variable in a parameter file. Instances of subclasses also 
   * contain a pointer to a C++ variable of the appropriate type.
   *
   * \ingroup Param_Module
   */
   class Parameter : public ParamComponent
   {

   public:

      /// Width of output field for a scalar variable.
      static const int Width = 20;

      /// Precision for io of floating point data field.
      static const int Precision = 12;

      /**
      * Constructor.
      *
      * \param label       label string preceding value in file format
      * \param isRequired  Is this a required parameter?
      */
      Parameter(const char *label, bool isRequired = true);

      /// Destructor.
      virtual ~Parameter();

      /// Return label string.
      std::string label() const;

      /// Is this an optional parameter?
      bool isRequired() const;

   protected:

      /// Label object that contains parameter label string.
      Label label_;

   };

}
#endif
