#ifndef TEXT_PARAM_COMPOSITE_H
#define TEXT_PARAM_COMPOSITE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>     // base class
#include <util/global.h>            

#include <vector>

namespace Util
{

   /**
   * A text representation of a ParamComposite.
   *
   * A TextComposite reads and writes a simple text representation of
   * a TextComposite file format, without interpreting its contents or
   * associating parameter values with variables. 
   *
   * \ingroup Param_Module
   */
   class TextComposite : public ParamComponent
   {
   
   public:
  
      /** 
      * Constructor
      */
      TextComposite();
  
      /** 
      * Virtual destructor.
      */
      virtual ~TextComposite();
   
      /**
      * Resets TextComposite to its empty state.
      */
      void resetParam();
   
      /** 
      * Read all parameters from an input stream.
      *
      * \param in input stream for reading
      */
      virtual void readParam(std::istream &in);
   
      /** 
      * Write all parameters to an output stream.
      *
      * This default implementation writes all parameters to file,
      * descending children recursively. 
      *
      * \param out output stream for reading
      */
      void writeParam(std::ostream &out);
   
   private:
   
      /// Array of pointers to ParamComponent objects.
      std::vector< std::string > lines_;    

   };
 
}
#endif
