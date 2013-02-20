#ifndef UTIL_BEGIN_H
#define UTIL_BEGIN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>

#include <string>

namespace Util
{

   /**
   * Beginning line of a composite parameter block.
   *
   * \ingroup Param_Module
   */
   class Begin : public ParamComponent
   {

   public:

      /**
      * Constructor.
      */
      Begin(const char* label);

      // Default destructor.

      /**
      * Read a comment
      *
      * \param in input stream
      */
      virtual void readParam(std::istream &in);

      /**
      * Read a comment
      *
      * \param out output stream
      */
      virtual void writeParam(std::ostream &out);

      /**
      * Return comment string
      */
      std::string string();

      /**
      * Do-nothing implementation of virtual resetParam function.
      */
      virtual void resetParam();

   private:

      /// Classname string
      std::string label_;

   };

   /*
   * Return comment string
   */
   inline std::string Begin::string()
   {  return label_; }

} 
#endif
