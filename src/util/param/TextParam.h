#ifndef TEXT_PARAM_H
#define TEXT_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>     // base class
#include <util/global.h>            

#include <vector>

namespace Util
{

   /**
   * A text representation of a Param.
   *
   * A TextParam reads and write a simple text representation of a parameter,
   * without interpreting its contents, or associating parameter values with
   * variables. 
   *
   * \ingroup Param_Module
   */
   class TextParam 
   {
   
   public:
  
      /** 
      * Constructor
      */
      TextParam();
  
      /** 
      * Virtual destructor.
      */
      virtual ~TextParam();
   
      /** 
      * Read all parameters from an input stream.
      *
      * \param in input stream for reading
      */
      virtual void read(&std::string* ptr);
   
      /*
      * Write all parameters to an output stream.
      *
      * This default implementation writes all parameters to file,
      * descending children recursively. 
      *
      * \param out output stream for reading
      */
      void write(std::ostream &out);
   
      //@}
      
   private:
  
      Label label_;

      std::string value_; 

   };
 
}
#endif
