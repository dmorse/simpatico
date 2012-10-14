#ifndef UTIL_END_H
#define UTIL_END_H

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
   * End bracket of a ParamComposite parameter block.
   *
   * \ingroup Param_Module
   */
   class End : public ParamComponent
   {

   public:

      /// Constructor.
      End();

      /// Destructor.
      virtual ~End();

      /**
      * Read the closing bracket.
      *
      * \param in input stream
      */
      virtual void readParam(std::istream &in);

      /**
      * Write the closing bracket.
      *
      * \param out output stream
      */
      virtual void writeParam(std::ostream &out);

      /// Do-nothing implementation of virtual resetParam function.
      void resetParam();

   };

} 
#endif
